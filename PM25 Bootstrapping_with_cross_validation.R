library(reshape2)
library(plyr)
library(mblm)
library(boot)
library(ggplot2)
library(gridExtra)
library(DAAG)


setwd("M:/AIR/MTSS/Python/Programs/PM25 Bootstrapping/Atlanta/Bootstrapping Results csa 122 2015-02-13")

seed = 12345

set.seed(seed) # Must set seed for reproducable results using bootstrapping

start_time = Sys.time()
print('------------- PM2.5 Incomplete Data Modeling and Bootstrapping ----------------')


# Define filenames
raw_fname = 'PM25 Daily Site Values_csa 122_2015-02-13.txt'
calendar_fname = 'aqs_sample_calendar.txt'
coll_freq_fname = 'site_collection_frequencies_2015-02-13.txt'
first_date = as.Date('2010-01-01')
last_date = as.Date('2014-12-31')
city_string = 'Atlanta'

### Run Parameters ###
min_sample_pairs = 100
min_quarterly_pairs = 20
sites_to_exclude = c()#c('13-067-0004', '13-089-2001', '13-121-0032') # Manually specify any sites to exclude from modeling analysis
bootstrap_runs = 1000
model_type = 'theil' #'theil' or 'least squares'
pm_naaqs_level = 12
num_folds = 10  # Number of folds to be used in k-folds evaluation of linear models
#i = 1 # testing code
######################

# Function Definitions

create_sample_calendar = function(site_freq, aqs_calendars, first_date, last_date){
    
  every_day_sample_days = aqs_calendar$aqs_date
  every_3rd_sample_days = aqs_calendar[aqs_calendar$nams_every_3rd_day_ind == 'Yes',
                                       ]$aqs_date
  every_6th_sample_days = aqs_calendar[aqs_calendar$nams_every_6th_day_ind == 'Yes',
                                       ]$aqs_date
  site_codes = unique(site_freq$site_code)
  # Set missing end dates to last_date
  site_freq$pri_mon_end_date[is.na(site_freq$pri_mon_end_date)] = last_date
  site_freq$req_coll_freq_end_date[is.na(site_freq$req_coll_freq_end_date)] = last_date
  # Calculate min end date and max begin date
  site_freq['min_end_date'] = as.Date(apply(site_freq[c('pri_mon_end_date','req_coll_freq_end_date')], 
                                            1, min, na.rm=T))
  site_freq['max_begin_date'] = as.Date(apply(site_freq[c('pri_mon_begin_date','req_coll_freq_begin_date')], 
                                              1, max, na.rm=T))
  # Set max begin date to first date if it is before first date
  site_freq$max_begin_date[site_freq$max_begin_date < first_date] = first_date
  
  scheduled_sample_days = NA # Initialize variable
  
  for(site_code in site_codes){
    coll_frequencies = site_freq[site_freq$site_code == site_code,]
    coll_frequencies = coll_frequencies[order(coll_frequencies$max_begin_date,
                                              decreasing=F),]
    # Build calendar for each site
    for(i in seq(nrow(coll_frequencies))){
      begin_date = coll_frequencies[i,]$max_begin_date
      end_date = coll_frequencies[i,]$min_end_date
      freq = coll_frequencies[i,]$cf_coll_freq_code
      if(freq == '1'){ # Daily samples
        dates = aqs_calendar[aqs_calendar$aqs_date >= begin_date &
                               aqs_calendar$aqs_date <= end_date,]$aqs_date
      }
      else if(freq == '3'){ # 1 in 3 day sampling
        dates = aqs_calendar[ aqs_calendar$nams_every_3rd_day_ind == 'Yes' &
                                (aqs_calendar$aqs_date >= begin_date &
                                   aqs_calendar$aqs_date <= end_date),]$aqs_date
      }
      else if(freq == '6'){ # 1 in 6 day sampling
        dates = aqs_calendar[ aqs_calendar$nams_every_6th_day_ind == 'Yes' &
                                (aqs_calendar$aqs_date >= begin_date &
                                   aqs_calendar$aqs_date <= end_date),]$aqs_date
      }
      else{
        stop(paste0('Error: unrecognized sampling frequency: ', freq))
      }
      # Add sample days to scheduled_sample_days dataframe
      if(!is.data.frame(scheduled_sample_days)){ # Initialize dataframe
        scheduled_sample_days = data.frame('site_code' = site_code,
                                           'sample_date' = dates)
      }
      else if (length(dates) > 0){
        #print(paste0('adding site: ', site_code)) #testing code
        rows_to_add = data.frame('site_code' = site_code,
                                 'sample_date' = dates)
        scheduled_sample_days = rbind(scheduled_sample_days, rows_to_add)
      }
      else{} # Pass
    }
  }
  return(scheduled_sample_days)
}

get_calendar_quarter = function(row){
  # Takes a vector of (sample_date, scheduled_date) and returns
  # the calendar quarter
  if(!(is.na(row[2]))){ date = row[2]} # Prefer scheduled date
  else { date = row[1]} # If scheduled date is NA, use sample date
  date = as.POSIXlt(date)
  month_number = as.numeric(format(date, '%m'))
  if(month_number %in% c(1,2,3)){quarter = 1}
  else if(month_number %in% c(4,5,6)){quarter = 2}
  else if(month_number %in% c(7,8,9)){quarter = 3}
  else if(month_number %in% c(10,11,12)){quarter = 4}
  else {stop(paste0('Month out of bounds: ', month_number))}  
}

compare_predictor_sites = function(pm_data, incomplete_sites, candidate_sites, seed=seed){
  # Takes a dataframe of daily pm_data and a list of incomplete sites and candidate sites,
  # and returns a dataframe of regression related statistics for each site pair of incomplete
  # site and candidate site
  
  # Initialize variables
  site_regressions = NA
  lm_residuals_df = NA
  theil_residuals_df = NA
  
  for(incomplete_site in incomplete_sites){
    for(candidate_site in candidate_sites){
      sample_pairs = merge( pm_data[pm_data$site_code == incomplete_site,][c('site_code', 'sample_date',
                                                                             'daily_value', 'quarter')],
                            pm_data[pm_data$site_code == candidate_site,][c('site_code', 'sample_date',
                                                                            'daily_value', 'quarter')],
                            by = c('sample_date', 'quarter'))
      sample_pairs = rename(sample_pairs, c('site_code.x' = 'incomplete_site',
                                            'site_code.y' = 'candidate_site',
                                            'daily_value.x' = 'incomplete_site_value',
                                            'daily_value.y' = 'candidate_site_value'))
      sample_pairs = sample_pairs[ complete.cases(sample_pairs),]
      n_pairs = nrow(sample_pairs)
      n_pairs_q1 = nrow(sample_pairs[sample_pairs$quarter == 1,])
      n_pairs_q2 = nrow(sample_pairs[sample_pairs$quarter == 2,])
      n_pairs_q3 = nrow(sample_pairs[sample_pairs$quarter == 3,])
      n_pairs_q4 = nrow(sample_pairs[sample_pairs$quarter == 4,])
      if(n_pairs >= min_sample_pairs & n_pairs_q1 >= min_quarterly_pairs & 
           n_pairs_q2 >= min_quarterly_pairs & n_pairs_q3 >= min_quarterly_pairs &
           n_pairs_q4 >= min_quarterly_pairs){
        min_pairs_met = 'Yes'
        
      }else{min_pairs_met = 'No'}
      
      lm_results = lm(data = sample_pairs, formula = incomplete_site_value ~ candidate_site_value)
      slope = lm_results$coefficients[2]
      intercept = lm_results$coefficients[1]
      r2 = summary(lm_results)$r.squared
      lm_residuals = lm_results$residuals
      
      theil_results = mblm(incomplete_site_value ~ candidate_site_value, dataframe = sample_pairs)
      theil_slope = theil_results$coefficients[2]
      theil_intercept = theil_results$coefficients[1]
      theil_residuals = theil_results$residuals
      
      pearson_r = cor(sample_pairs$incomplete_site_value, sample_pairs$candidate_site_value,
                      method='pearson')
      spearman_r = cor(sample_pairs$incomplete_site_value, sample_pairs$candidate_site_value,
                       method='spearman')
      
      #print(paste0('Incomplete Site: ', incomplete_site, ' Candidate Site: ', candidate_site))
      
      # Cross-validation of linear models using k-folds
      lm_cross_validation = CVlm(df = sample_pairs, form.lm=lm_results, m=num_folds, plot= T, printit=F)
      # Extract Mean Squared Predition Error
      lm_mse = attributes(lm_cross_validation)$ms
      # Calculate Root Mean Squared Error (RMSE)
      lm_rmse = sqrt(lm_mse)
      
      # Change class of mblm results from c('mblm', 'lm') to 'lm' so that CVlm function
      # won't throw an error.  mblm results are designed to be compatable with lm methods.
      class(theil_results) = 'lm' 
      
      theil_cross_validation = CVlm(df = sample_pairs, form.lm=theil_results, m=num_folds,
                                    plot=T, printit=F)
      # Extract Mean Squared Predition Error
      theil_mse = attributes(theil_cross_validation)$ms
      # Calculate Root Mean Squared Error (RMSE)
      theil_rmse = sqrt(theil_mse)
      
      if(!is.data.frame(site_regressions)){ # Create dataframe of model results
        site_regressions = data.frame('incomplete_site' = incomplete_site,
                                      'candidate_site' = candidate_site,
                                      'n_pairs' = n_pairs,
                                      'n_pairs_q1' = n_pairs_q1,
                                      'n_pairs_q2' = n_pairs_q2,
                                      'n_pairs_q3' = n_pairs_q3,
                                      'n_pairs_q4' = n_pairs_q4,
                                      'min_pairs_met' = min_pairs_met,
                                      'slope' = slope,
                                      'intercept' = intercept,
                                      'r2' = r2,
                                      'theil_slope' = theil_slope,
                                      'theil_intercept' = theil_intercept,
                                      'pearson_r' = pearson_r,
                                      'spearman_r' = spearman_r,
                                      'lm_rmse' = lm_rmse,
                                      'theil_rmse' = theil_rmse)
      }else{ # Append row to site_regressions
        row = data.frame('incomplete_site' = incomplete_site,
                         'candidate_site' = candidate_site,
                         'n_pairs' = n_pairs,
                         'n_pairs_q1' = n_pairs_q1,
                         'n_pairs_q2' = n_pairs_q2,
                         'n_pairs_q3' = n_pairs_q3,
                         'n_pairs_q4' = n_pairs_q4,
                         'min_pairs_met' = min_pairs_met,
                         'slope' = slope,
                         'intercept' = intercept,
                         'r2' = r2,
                         'theil_slope' = theil_slope,
                         'theil_intercept' = theil_intercept,
                         'pearson_r' = pearson_r,
                         'spearman_r' = spearman_r,
                         'lm_rmse' = lm_rmse,
                         'theil_rmse' = theil_rmse)
        site_regressions = rbind(site_regressions, row)
      }
      # Create or add to residuals dataframes
      if(!is.data.frame(lm_residuals_df)){
          lm_residuals_df = data.frame('incomplete_site' = incomplete_site,
                                       'candidate_site' = candidate_site,
                                       'residuals' = lm_residuals)
      }else{
          rows = data.frame('incomplete_site' = incomplete_site,
                            'candidate_site' = candidate_site,
                            'residuals' = lm_residuals)
          lm_residuals_df = rbind(lm_residuals_df, rows)
      }
        
      if(!is.data.frame(theil_residuals_df)){
        theil_residuals_df = data.frame('incomplete_site' = incomplete_site,
                                       'candidate_site' = candidate_site,
                                       'residuals' = theil_residuals)
      }else{
        rows = data.frame('incomplete_site' = incomplete_site,
                            'candidate_site' = candidate_site,
                            'residuals' = theil_residuals)
        theil_residuals_df = rbind(theil_residuals_df, rows)
      }
        
    } 
  }
  return(list('site_regressions' = site_regressions,
              'lm_residuals' = lm_residuals_df,
              'theil_residuals' = theil_residuals_df))
}

select_best_models = function(site_regressions, model_type){
  # Returns a table of the best fit model (highest r squared for least squares, or spearman r for theil)
  # for each incomplete site. Excludes any predictor sites with not enough sample pairs
  incomplete_sites = unique(site_regressions$incomplete_site)
  candidate_regressions = site_regressions[ site_regressions$min_pairs_met == 'Yes',]
  best_models = NA # Initialize variable
  for(incomplete_site in incomplete_sites){
    models = candidate_regressions[ candidate_regressions$incomplete_site == incomplete_site,]
    if(model_type == 'least_squares'){
      best_model = models[ models$r2 == max(models$r2),]
    }
    else if(model_type == 'theil'){
      best_model = models[ models$spearman_r == max(models$spearman_r),]
    }
    else{
      stop(paste0('Unrecognized model type: ', model_type))
    }
    
    if(!is.data.frame(best_models)){ # Initialize dataframe
      best_models = best_model
    }
    else{
      best_models = rbind(best_models, best_model)
    }
  }
  return(best_models)
}

create_regression_plots = function(best_models, pm_data, model_type){

  plots = list()
  max_value = 0
  #equations = list()
  for(i in seq(nrow(best_models))){
    model = best_models[i,]
    incomplete_site = as.character(model$incomplete_site)
    candidate_site = as.character(model$candidate_site)
    sample_pairs = merge( pm_data[pm_data$site_code == incomplete_site,][c('site_code', 'sample_date',
                                                                           'daily_value', 'quarter')],
                          pm_data[pm_data$site_code == candidate_site,][c('site_code', 'sample_date',
                                                                          'daily_value', 'quarter')],
                          by = c('sample_date', 'quarter'))
    sample_pairs = rename(sample_pairs, c('site_code.x' = 'incomplete_site',
                                          'site_code.y' = 'candidate_site',
                                          'daily_value.x' = 'incomplete_site_value',
                                          'daily_value.y' = 'candidate_site_value'))
    sample_pairs = sample_pairs[ complete.cases(sample_pairs),]
    plot_name = paste0('x_',gsub('-', '', candidate_site), '_y_',
                       gsub('-', '', incomplete_site))
    sample_pairs_max = max( max(sample_pairs$incomplete_site_value), 
                            max(sample_pairs$candidate_site_value))
    if( sample_pairs_max > max_value){ max_value = sample_pairs_max}

    if( model_type == 'theil'){
      slope = model$theil_slope
      intercept = model$theil_intercept
      correlation = model$spearman_r
      equation = as.character(as.expression(
        substitute(italic(y) == slope ~ italic(x) + intercept*", "~rho~"="~correlation,
                                          list(slope = round(slope,2),
                                               intercept = round(intercept,2),
                                               correlation = round(correlation,2)))))
      #print(equation) # testing code
      ymax = max(sample_pairs$incomplete_site_value)
      plot = ggplot(data = sample_pairs, aes(x=candidate_site_value, y=incomplete_site_value)) +
        geom_point(color='navy', alpha=0.8) +
        xlab(candidate_site) + ylab(incomplete_site) +
        geom_abline(slope = 1, intercept = 0, linetype = 'longdash', color = 'black') +
        geom_abline(slope = slope, intercept = intercept, color = 'red') +
        #geom_text(aes(x=5, y=27, label = eval(equation)),parse=T)
        annotate("text", x=5, y=ymax*1.05, label=equation, parse=T, hjust=0)
    
      }else if( model_type == 'least_squares'){
        slope = model$slope
        intercept = model$intercept
        correlation = model$r2
        equation = as.character(as.expression(
          substitute(italic(y) == slope ~ italic(x) + intercept*", "~italic(r^2)~"="~correlation,
                     list(slope = round(slope,2),
                          intercept = round(intercept,2),
                          correlation = round(correlation,2)))))
        
        plot = ggplot(data = sample_pairs, aes(x=candidate_site_value, y=incomplete_site_value)) +
          geom_point() +
          xlab(candidate_site) + ylab(incomplete_site) +
          geom_abline(slope = 1, intercept = 0, linetype = 'longdash', color = 'blue') +
          geom_abline(slope = slope, intercept = intercept, color = 'red') +
          geom_text(aes(x=5, y=27, label = equation),parse=T)
      }else{
        stop(paste0('Unrecognized model type: ', model_type))
      }
    plots[[plot_name]] = plot
  }
  plots[['axis_max']] = max_value * 1.05 # Use ggplot default additive axis factor of .05\
  
  return(plots)
}

create_seasonal_model_test_data = function(best_models, pm_data){
  
  model_test_data = NA
  for(i in seq(nrow(best_models))){
    model = best_models[i,]
    incomplete_site = as.character(model$incomplete_site)
    predictor_site = as.character(model$candidate_site)
    sample_pairs = merge( pm_data[pm_data$site_code == incomplete_site,][c('site_code', 'sample_date',
                                                                           'daily_value', 'quarter')],
                          pm_data[pm_data$site_code == predictor_site,][c('site_code', 'sample_date',
                                                                          'daily_value', 'quarter')],
                          by = c('sample_date', 'quarter'))
    sample_pairs = rename(sample_pairs, c('site_code.x' = 'incomplete_site',
                                          'site_code.y' = 'predictor_site',
                                          'daily_value.x' = 'incomplete_site_value',
                                          'daily_value.y' = 'predictor_site_value'))
    sample_pairs = sample_pairs[ complete.cases(sample_pairs),]
    
    lm_results = lm(data = sample_pairs, formula = incomplete_site_value ~ predictor_site_value)
    lm_slope = lm_results$coefficients[2]
    lm_intercept = lm_results$coefficients[1]
    r2 = summary(lm_results)$r.squared
    lm_residuals = lm_results$residuals
    
    theil_results = mblm(incomplete_site_value ~ predictor_site_value, dataframe = sample_pairs)
    theil_slope = theil_results$coefficients[2]
    theil_intercept = theil_results$coefficients[1]
    theil_residuals = theil_results$residuals
    
    ### Analyze model error by quarter
    # This analysis is included to investigate model performance on time series data,
    # which is not actually iid as linear modeling, bootstrapping, and k-folds cross validation assume.
    
    model_test = sample_pairs
    model_test['lm_modeled_value'] = model_test$predictor_site_value * lm_slope + lm_intercept
    model_test['lm_residual'] = model_test$incomplete_site_value - model_test$lm_modeled_value
    model_test['theil_modeled_value'] = model_test$predictor_site_value * theil_slope + theil_intercept
    model_test['theil_residual'] = model_test$incomplete_site_value - model_test$theil_modeled_value 
    
    model_test['year'] = as.POSIXlt(model_test$sample_date)$year + 1900
    model_test['year_quarter'] = paste0(model_test$year, '-Q', model_test$quarter)
    
    if(is.na(model_test_data)){
      model_test_data = model_test
    }else{
      model_test_data = rbind(model_test_data, model_test)
    }
  }
  return(model_test_data)
}


model_missing_data = function(data, data_comp, best_models, model_type){
  # Model missing data in place in data frame, and return data frame with new column of modeled values
  # Only runs model for missing values at incomplete sites during incomplete quarters
  
  incomplete_sites = best_models$incomplete_site
  # Create placeholder columns for modeled values and data source
  data['predictor_site'] = NA
  data['predictor_site_value'] = NA
  data['modeled_value'] = NA 
  data['data_source'] = 'sample' # default data source is from sample record
  
  # Run model in place in recent_data
  for(i in seq(nrow(data))){
    row = data[i,]
    site = as.character(row$site_code)
    year_quarter = row$year_quarter
    daily_value = row$daily_value
    sample_date = row$sample_date
    if(site %in% incomplete_sites & is.na(daily_value)){ # Check for missing value at incomplete site
      incomplete_quarters = data_comp[ data_comp$site_code == site &
                                         data_comp$year %in% most_recent_dv_period &
                                         (data_comp$pct_complete < 75 | 
                                            is.na(data_comp$pct_complete)),]$year_quarter
      if(year_quarter %in% incomplete_quarters){ # Check if this missing value is in an incomplete quarter
        # Perform modeling of the incomplete value
        model_results = best_models[ best_models$incomplete_site == site,]
        predictor_site = as.character(model_results$candidate_site)
        # Select slope and intercept to use based on model_type
        if(model_type == 'least_squares'){
          slope = model_results$slope
          intercept = model_results$intercept
        }
        else if(model_type == 'theil'){
          slope = model_results$theil_slope
          intercept = model_results$theil_intercept
        }
        else{
          stop(paste0('Unrecognized model type: ', model_type))
        }
        # Find data value from predictor site
        predictor_site_value = data[ data$site_code == predictor_site &
                                       data$sample_date == sample_date,]$daily_value
        data[i,]$predictor_site = predictor_site
        data[i,]$predictor_site_value = predictor_site_value
        
        # Add modeled value into data
        modeled_value = (slope * predictor_site_value) + intercept
        data[i,]$modeled_value = modeled_value
        data[i,]$data_source = 'model'
      }
    }
    else{ # Data from sample record kept and go to next row in recent_data
      data[i,]$data_source = 'sample'
    } 
  }
  # Add column of combined dataset with collected and modeled values
  data['combined'] = data$daily_value
  data[is.na(data$daily_value),]$combined = 
    data[is.na(data$daily_value),]$modeled_value
  #data[is.na(data$combined),]$data_source = NA
  return(data)
}

calc_annual_dv = function(data){
  # Takes in a data table of three years of data from a site, and returns a one row
  # data frame with the site code, dv_year, complete quarters, and annual dv
  #data = incomplete_site_data
  dv_year = max(data$year)
  site_code = unique(data$site_code)
  
  # Check data completeness
  completeness = ddply(data, ~site_code + year + quarter + year_quarter, summarize,
                       scheduled_samples = length(unique(na.exclude(scheduled_date))),
                       actual_values = length(na.exclude(na.exclude(value))),
                       pct_complete = round((actual_values/scheduled_samples)*100,1))
  complete_quarters = nrow(completeness[ completeness$pct_complete >= 75,])
  
  if(complete_quarters == 12){ # Calculate DV if all quarters are complete
    # Calculate quarterly means, annual means, and 3-year DV
    quarterly_means = ddply(data, ~site_code + year + quarter + year_quarter, summarize,
                            quarterly_mean = mean(value, na.rm=T))
    annual_means = ddply(quarterly_means, ~site_code + year, summarize,
                         annual_mean = mean(quarterly_mean, na.rm=T))
    annual_dv = round(mean(annual_means$annual_mean, na.rm=T), 2)
  }
  else{  annual_dv = NA }
  
  result = data.frame('site_code' = site_code,
                      'dv_year' = dv_year,
                      'complete_quarters' = complete_quarters,
                      'annual_dv' = annual_dv)
  return(result)
  
}

dv_for_bootstrap = function(residuals, indices, site_data){
  # Function to provide to boot function as 'statistic' argument
  # Boot function will randomly select residuals to apply to modeled values
  # in pm_data using indices, and then run calc_annual_dv to generate a design_value
  
  random_residuals = residuals[indices] # allows boot to select residuals to apply
  # Select only as many residuals to apply as rows in pm_data
  residuals_to_apply = 
    random_residuals[1:nrow(site_data[site_data$data_source == 'model' &
                                        !is.na(site_data$data_source),])]
  #print('Residuals to apply:') # testing code
  #print(head(residuals_to_apply)) #testing code
  #print('Data before residuals applied: ') # testing code
  #print(head(site_data[site_data$data_source == 'model' & # testing code
  #                      !is.na(site_data$data_source),]$value))
  
  # Apply random residuals only to modeled values in dataset
  new_data = site_data
  new_data[new_data$data_source == 'model' & 
             !is.na(new_data$data_source),]$value = 
    new_data[new_data$data_source == 'model' & 
               !is.na(new_data$data_source),]$value + residuals_to_apply
  # Calculate annual design value using new_data
  
  #print('Data after residuals applied: ') # testing code
  #print(head(new_data[new_data$data_source == 'model' & # testing code
  #                      !is.na(new_data$data_source),]$value))
  
  annual_dv = calc_annual_dv(data=new_data)$annual_dv
  #print(paste0('DV: ', annual_dv))
  #print('------------------------------------------------------')
  return(annual_dv)
}


# Import Data
print('Importing data...')
pm_raw = read.table(raw_fname, header=T, sep='\t', quote="", na.strings='None',
                    colClasses='character')
pm_raw$daily_value = as.numeric(pm_raw$daily_value)
pm_raw$sample_day = as.Date(pm_raw$sample_day, format='%Y-%m-%d %H:%M:%S')
pm_raw$scheduled_date = as.Date(pm_raw$scheduled_date, format='%Y-%m-%d %H:%M:%S')

aqs_calendar = read.table(calendar_fname, header=T, sep='\t', quote="",
                          na.strings='None', colClasses='character')
aqs_calendar$aqs_date = as.Date(aqs_calendar$aqs_date, format='%Y-%m-%d %H:%M:%S')

site_coll_freq = read.table(coll_freq_fname, header=T, sep='\t', quote="",
                            na.strings='None', colClasses='character')
site_coll_freq$pri_mon_begin_date = as.Date(site_coll_freq$pri_mon_begin_date,
                                            format = '%Y-%m-%d %H:%M:%S')
site_coll_freq$pri_mon_end_date = as.Date(site_coll_freq$pri_mon_end_date,
                                            format = '%Y-%m-%d %H:%M:%S')
site_coll_freq$req_coll_freq_begin_date = as.Date(site_coll_freq$req_coll_freq_begin_date,
                                            format = '%Y-%m-%d %H:%M:%S')
site_coll_freq$req_coll_freq_end_date = as.Date(site_coll_freq$req_coll_freq_end_date,
                                            format = '%Y-%m-%d %H:%M:%S')
# Remove manually selected sites to exclude
site_coll_freq = site_coll_freq[ !(site_coll_freq$site_code %in% sites_to_exclude),]

# Create table of scheduled sample days by site
scheduled_sample_days = create_sample_calendar(site_coll_freq, aqs_calendars,
                                               first_date, last_date)

# Merge PM raw data with scheduled sample days

pm_data  = merge(x=scheduled_sample_days, 
                 y=pm_raw[c('site_code', 'sample_day', 'daily_value', 
                            'scheduled_date', 'creditable_day_indicator')],
                 by.x = c('site_code', 'sample_date'),
                 by.y = c('site_code', 'sample_day'),
                 all= T)
# To correctly record scheduled dates for missed samples, copy sample_date into
# scheduled date for missed samples (records where scheduled date and value are NA)
pm_data[is.na(pm_data$scheduled_date) & is.na(pm_data$daily_value),]$scheduled_date = 
  pm_data[is.na(pm_data$scheduled_date) & is.na(pm_data$daily_value),]$sample_date


# Use scheduled date, not sample date to correctly calculate 
# calendar quarters for makeup samples
pm_data['quarter'] = apply(pm_data, 1, get_calendar_quarter) 

pm_data['year'] = NA # Initialize column
pm_data[!is.na(pm_data$scheduled_date),]$year =
  1900 + as.POSIXlt(pm_data[!is.na(pm_data$scheduled_date),]$scheduled_date)$year
pm_data[is.na(pm_data$scheduled_date),]$year =
  1900 + as.POSIXlt(pm_data[is.na(pm_data$scheduled_date),]$sample_date)$year

pm_data['year_quarter'] = paste0(pm_data$year, '-Q', pm_data$quarter)
# Exclude data from manually selected sites
pm_data = pm_data[ !(pm_data$site_code %in% sites_to_exclude),]



# Calculate data completeness
print('Calculating data completeness...')
pm_data_comp = ddply(pm_data, ~site_code + year + quarter + year_quarter, summarize,
                     scheduled_samples = length(unique(na.exclude(scheduled_date))),
                     creditable_samples = length(na.exclude(creditable_day_indicator)),
                     pct_complete = round((creditable_samples/scheduled_samples)*100,1))
#crosstab_data_comp = dcast(pm_data_comp[c('site_code', 'year_quarter', 'pct_complete')],
#                          year_quarter ~ site_code, value.var = 'pct_complete')
crosstab_data_comp = dcast(pm_data_comp[c('site_code', 'year_quarter', 'pct_complete')],
                           site_code ~ year_quarter, value.var = 'pct_complete')
data_comp = melt(crosstab_data_comp, id.vars=c('site_code'))
data_comp = rename(data_comp, c('variable'='year_quarter',
                                'value' = 'pct_complete'))
data_comp['year'] = as.numeric(substr(data_comp$year_quarter, 1,4))

#write.csv(pm_data_comp, 'pm_data_comp.csv', row.names=F)
write.csv(crosstab_data_comp, 'crosstab_data_comp.csv', row.names=F)

# Identify incomplete sites and candidate sites for modeling incomplete data
last_year = 1900 + as.POSIXlt(last_date)$year
most_recent_dv_period = c(last_year, last_year-1, last_year-2)

complete_quarter_count = ddply(pm_data_comp[pm_data_comp$year %in% most_recent_dv_period,], 
                               ~site_code, summarize,
                               complete_quarter_count = sum(pct_complete >= 75))
complete_sites = unique(complete_quarter_count[ 
  complete_quarter_count$complete_quarter_count == 12,]$site_code)
sites = unique(pm_data_comp$site_code)
incomplete_sites = as.character(sites[!sites %in% complete_sites])


##### Modeling of incomplete data ######
# For purposes of Louisville analysis, only use sites with complete design values
# for the most recent design value period
print('Running regressions between incomplete sites and all candidate sites...')
candidate_sites = complete_sites
compare_predictor_site_results = 
  compare_predictor_sites(pm_data, incomplete_sites, candidate_sites)
site_regressions = compare_predictor_site_results$site_regressions
lm_residuals = compare_predictor_site_results$lm_residuals
theil_residuals = compare_predictor_site_results$theil_residuals
write.csv(site_regressions, 'site_regressions.csv', row.names=F)
print('Regressions calculated and written to site_regressions.csv')

print('Selecting best-fit models for each incomplete site')
best_models = select_best_models(site_regressions, model_type)
write.csv(best_models, 'best_fit_models.csv', row.names=F)
print('Best-fit models selected and written to best_fit_models.csv')

# Create regression plots of best-fit models
regression_plots = create_regression_plots( best_models, pm_data, model_type)
print('Regression plots created.')

# Select only pm_data and data_comp in most_recent_dv_period for use in model
print('Modeling incomplete data using best-fit models...')
recent_data =  pm_data[ pm_data$year %in% most_recent_dv_period,]
recent_data_comp = data_comp[ data_comp$year %in% most_recent_dv_period,]

results = model_missing_data(data = recent_data,
                             data_comp = recent_data_comp,
                             best_models = best_models,
                             model_type = model_type)

print(paste0('Incomplete data modeling finished. Total values modeled: ',
             length(na.exclude(results$modeled_value))))

# Calculate new design values for incomplete sites based on the combination of
# sampled and modeled data
print("Calculating new design values based on combined sampled and modeled data...")
modeled_dvs = NA # Initialize variable
for(incomplete_site in incomplete_sites){
  # Select data for submittal to design value function
  incomplete_site_data = results[ results$site_code == incomplete_site,]
  incomplete_site_data = incomplete_site_data[c('site_code', 'sample_date', 
                                                'scheduled_date', 'combined',
                                                'quarter', 'year', 'year_quarter')]
  incomplete_site_data = rename(incomplete_site_data, c('combined' = 'value'))
  if(!is.data.frame(modeled_dvs)){
    modeled_dvs = calc_annual_dv(incomplete_site_data)
  }
  else{
    row = calc_annual_dv(incomplete_site_data)
    modeled_dvs = rbind(modeled_dvs, row)
  }
}
write.csv(modeled_dvs, 'modeled_dvs.csv', row.names=F)
print('Modeled design values written to modeled_dvs.csv')

# Perform bootstrapping to estimate variance of modeled design values
print('Performing bootstrapping to estimate variance of modeled design values...')
if(model_type == 'least_squares'){
  residuals = lm_residuals
}else if(model_type == 'theil'){
  residuals = theil_residuals
}else{
  stop(paste0('Unrecognized model type: ', model_type))
}

boot_results_list = list()
for(incomplete_site in incomplete_sites){
  candidate_site = as.character(best_models[ best_models$incomplete_site == 
                                  incomplete_site,]$candidate_site)
  site_residuals = residuals[ residuals$incomplete_site == incomplete_site &
                                residuals$candidate_site == candidate_site,]$residuals
  incomplete_site_data = results[ results$site_code == incomplete_site,]
  incomplete_site_data = incomplete_site_data[c('site_code', 'sample_date', 
                                                'scheduled_date', 'combined', 'data_source',
                                                'quarter', 'year', 'year_quarter')]
  incomplete_site_data = rename(incomplete_site_data, c('combined' = 'value'))
  boot_results = boot(data=site_residuals, statistic=dv_for_bootstrap, 
                      R=bootstrap_runs, site_data=incomplete_site_data)
  boot_results_list[[incomplete_site]] = boot_results
}

# Summarize Bootstrap Results
print('Summarizing bootstrap results...')
bootstrapped_dvs = NA
bootstrap_summary = NA
for(i in seq(length(boot_results_list))){
  # Combine bootstrapped dv results into one table for graphing
  incomplete_site = names(boot_results_list[i])
  b_result = boot_results_list[[i]]
  if(!is.data.frame(bootstrapped_dvs)){
    bootstrapped_dvs = data.frame('incomplete_site' = incomplete_site,
                                  'design_value'= b_result$t)
  }else{
    rows = data.frame('incomplete_site' = incomplete_site,
                                         'design_value'= b_result$t)
    bootstrapped_dvs = rbind(bootstrapped_dvs, rows)
  }
  # Summarize results of bootstrapping
  predictor_site = best_models[ best_models$incomplete_site == 
                                  incomplete_site,]$candidate_site
  modeled_dv = modeled_dvs[ modeled_dvs$site_code == incomplete_site,]$annual_dv
  t0 = b_result$t0
  if(is.na(t0)){
    conf_interval_lower_bound_95 = NA
    conf_interval_upper_bound_95 = NA
    max_bootstrapped_dv = NA
    min_bootstrapped_dv = NA
  }else{
    boot_ci = boot.ci(b_result, type = 'basic')
    conf_interval_lower_bound_95 = boot_ci$basic[4]
    conf_interval_upper_bound_95 = boot_ci$basic[5]
    max_bootstrapped_dv = max(b_result$t)
    min_bootstrapped_dv = min(b_result$t)
  }
  
  if(!is.data.frame(bootstrap_summary)){
    bootstrap_summary = data.frame('incomplete_site' = incomplete_site,
                                   'predictor_site' = predictor_site,
                                   'bootstrap_runs' = length(b_result$t),
                                   'modeled_dv' = modeled_dv,
                                   'dv_t0' = t0,
                                   'conf_interval_lower_bound_95' = conf_interval_lower_bound_95,
                                   'conf_interval_upper_bound_95' = conf_interval_upper_bound_95,
                                   'max_bootstrapped_dv' = max_bootstrapped_dv,
                                   'min_bootstrapped_dv' = min_bootstrapped_dv)
  }else{
    row = data.frame('incomplete_site' = incomplete_site,
                     'predictor_site' = predictor_site,
                     'bootstrap_runs' = length(b_result$t),
                     'modeled_dv' = modeled_dv,
                     'dv_t0' = t0,
                     'conf_interval_lower_bound_95' = conf_interval_lower_bound_95,
                     'conf_interval_upper_bound_95' = conf_interval_upper_bound_95,
                     'max_bootstrapped_dv' = max_bootstrapped_dv,
                     'min_bootstrapped_dv' = min_bootstrapped_dv)
    bootstrap_summary = rbind(bootstrap_summary, row)
  }
}

write.csv(bootstrap_summary, 'bootstrap_summary.csv')

bootstrapped_dvs = bootstrapped_dvs[complete.cases(bootstrapped_dvs),] # exclude NAs

histogram_plot = ggplot(data=bootstrapped_dvs, aes(x=design_value, y=..density..)) + # 
  geom_bar(binwidth=0.01, fill='cornflowerblue', position='identity') +
  geom_density() + 
  theme_bw() +
  facet_grid(incomplete_site~.) +
  ggtitle(paste0(city_string, ' Bootstrap Results of Modeled 2011-2013 PM2.5 Annual Design Values')) +
  xlab(as.expression('Modeled Design Value ('~mu~'g/'~m^3~')')) +
  ylab(paste0('Frequency of Outcome in ', bootstrap_runs, ' Bootstrap Runs')) +
  geom_vline(xintercept = pm_naaqs_level, color = 'red', linetype = 'longdash') +
  scale_x_continuous(breaks = seq(9.0, 12.0, 0.2))

data_comp_heatmap = ggplot(data=data_comp, 
                           aes(x=year_quarter, y=site_code, fill=pct_complete)) +
  geom_tile() + 
  scale_fill_gradient2(name = "Percent Complete", low="orangered4", mid="white", 
                       high="dodgerblue4", midpoint = 75) +
  theme_minimal() + 
  ggtitle(paste0(city_string, ' Quarterly PM2.5 Data Completeness 2009-2013')) +
  xlab('Year and Quarter') + ylab('Site Code') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(city_string, ' Bootstrap Results 09-13.png'), histogram_plot,
       height = 7.5, width = 10, units='in')
ggsave(paste0(city_string, ' Data Completeness 09-13.png'), data_comp_heatmap,
       height = 7.5, width = 10, units='in')

### Seasonal Cross Validation of Linear models
model_test_data = create_seasonal_model_test_data(best_models, pm_data)
model_test_data ['site_pair'] = paste0(model_test_data$incomplete_site, ' / ', model_test_data$predictor_site)
  
# Evaluate model error by quarter
# Equivalent to a k-folds analysis, but using quarters instead of randomly assigned folds
# since the data to be modeled is often in a block rather than randomly distributed
if(model_type == 'lm'){
  model_error_by_quarter = ddply(model_test_data, ~incomplete_site + predictor_site + 
                                            site_pair + year + quarter + year_quarter, summarize,
                                          n = length(na.exclude(incomplete_site_value)),
                                          rmse = sqrt(mean(lm_residual^2)),
                                          median_residual = median(lm_residual))
}else if(model_type == 'theil'){
  model_error_by_quarter = ddply(model_test_data, ~incomplete_site + predictor_site + 
                                            site_pair + year + quarter + year_quarter, summarize,
                                          n = length(na.exclude(incomplete_site_value)),
                                          rmse = sqrt(mean(theil_residual^2)),
                                          median_residual = median(theil_residual))
}

model_error_timeseries = ggplot(data=model_error_by_quarter, 
                                aes(x=year_quarter, y=rmse,
                                    color=site_pair)) +
  geom_point() + geom_line(aes(group=site_pair)) +
  ylab('Root Mean Square Error of Model (ug/m3)') +
  xlab('Year and Quarter') + 
  ggtitle(paste0(city_string, 'Model Prediction Error by Calendar Quarter'))

median_residual_timeseries = ggplot(data=model_error_by_quarter, 
                                    aes(x=year_quarter, y=median_residual,
                                        color=site_pair)) +
  geom_point() + geom_line(aes(group=site_pair)) +
  ylab('Median Model Residual Value (ug/m3)') +
  xlab('Year and Quarter') + 
  ggtitle(paste0(city_string, 'Model Median Residual Value by Calendar Quarter'))

quarter_labeller = function(variable, value){
  quarter_names = c('Quarter 1' = 1, 'Quarter 2' = 2,
                    'Quarter 3' = 3, 'Quarter 4' = 4)
  return(names(quarter_names[quarter_names == value]))
}

model_error_bar = ggplot(data=model_error_by_quarter, aes(x=year, y=rmse, fill=site_pair)) +
  geom_bar(stat='identity', position='dodge') + 
  facet_grid(~quarter, labeller=quarter_labeller) +
  ylab('Root Mean Square Error of Model (ug/m3)') +
  xlab('Year') + 
  ggtitle(paste0(city_string, ' Model Prediction Error by Calendar Quarter'))+
  scale_fill_discrete(name="Incomplete Site/\nPredictor Site")

median_residual_bar = ggplot(data=model_error_by_quarter, aes(x=year, y=median_residual, fill=site_pair)) +
  geom_bar(stat='identity', position='dodge') + 
  facet_grid(~quarter) +
  ylab('Median Model Residual Value (ug/m3)') +
  xlab('Year') + 
  ggtitle(paste0(city_string, 'Model Median Residual Value by Calendar Quarter')) 

# Save Regression Plots
# Note: this code must be updated for each area to correctly save the plots
# using grid.arange

# Create plots with adjusted axes
regression_plots_adjusted = list()
axis_max = regression_plots$axis_max
for( i in seq(length(regression_plots)-1)){
  myplot = regression_plots[[i]]
  myname = names(regression_plots)[i]
  #print(myname)
  myplot = myplot + scale_x_continuous(limits = c(-5, axis_max)) +
    scale_y_continuous(limits = c(-5, axis_max))
  regression_plots_adjusted[[myname]] = myplot
}




end_time = Sys.time()
elapsed_time = end_time - start_time
print( paste0('Done! Elapsed Time: ', elapsed_time, ' minutes'))

# Export regression plots.  This code must be updated with each run!

png(paste0(city_string,'_Regressions.png'), width=20, height=15, units='in', res=150)
grid.arrange(regression_plots_adjusted$x_450370001_y_132450091,
             ncol=1, 
             main=paste0(city_string, ' PM2.5 Theil Regressions and Spearman Correlations 2009-2013'))
dev.off()


###################################################
# Run models for each incomplete site sequentially
#for( i in seq(nrow(best_models))){
#  model_results = best_models[i,]
#  incomplete_site = as.character(model_results$incomplete_site)
#  predictor_site = as.character(model_results$candidate_site)
  # Determine which calendar quarters are incomplete
  # Select only the incomplete and predictor site data for quarters that are 
  # incomplete at the incomplete site. Quarters with >75% data are not modeled
#  incomplete_quarters = data_comp[ data_comp$site_code == incomplete_site &
#                                     data_comp$year %in% most_recent_dv_period &
#                                     (data_comp$pct_complete < 75 | 
#                                     is.na(data_comp$pct_complete)),]$year_quarter
#  incomplete_and_predictor_site_data = pm_data[ pm_data$site_code %in% c(incomplete_site, predictor_site) &
#                                    pm_data$year %in% most_recent_dv_period &
#                                      pm_data$year_quarter %in% incomplete_quarters,]
#  incomplete_and_predictor_site_data['site_type'] = NA
#  if( incomplete_site %in% unique(incomplete_and_predictor_site_data$site_code)){ #check to see if any data exists for incomplete site
#    incomplete_and_predictor_site_data[ incomplete_and_predictor_site_data$site_code ==
#                                        incomplete_site,]$site_type = 'incomplete_site'}
#  incomplete_and_predictor_site_data[ incomplete_and_predictor_site_data$site_code ==
#                                        predictor_site,]$site_type = 'predictor_site'
#  input_data = dcast(data = incomplete_and_predictor_site_data,
#                     formula = sample_date ~ site_type,
#                     value.var = 'daily_value')
  # If there is no data from the incomplete site, add a column of NAs
#  if(!('incomplete_site' %in% colnames(input_data))){
#    input_data['incomplete_site'] = NA
#  }
#}

#> a = data.frame('letter' = c('a','b','c'), 'number' = c(1,2,3))
#> a$number = a$number + c(0.1,0.2,0.3)
######
#incomplete_site = incomplete_sites[5]
#candidate_site = as.character(best_models[ best_models$incomplete_site == 
#                                             incomplete_site,]$candidate_site)
#site_residuals = residuals[ residuals$incomplete_site == incomplete_site &
#                              residuals$candidate_site == candidate_site,]$residuals
#incomplete_site_data = results[ results$site_code == incomplete_site,]
#incomplete_site_data = incomplete_site_data[c('site_code', 'sample_date', 
#                                              'scheduled_date', 'combined', 'data_source',
#                                              'quarter', 'year', 'year_quarter')]
#incomplete_site_data = rename(incomplete_site_data, c('combined' = 'value'))
#boot_results = boot(data=site_residuals, statistic=dv_for_bootstrap, 
#                    R=bootstrap_runs, site_data=incomplete_site_data)

