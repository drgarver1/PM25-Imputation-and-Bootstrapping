# pm25_bootstrap.py
# This program runs the imputation and bootstrapping protocol on the incomplete
# PM2.5 sites in a CSA

from sqlalchemy import *
#from sqlalchemy.ext.sqlsoup import SqlSoup
import sqlsoup
import sys
import os
import datetime
sys.path.append('M:\AIR\MTSS\Python\Modules')
import mtss

starttime = datetime.datetime.now()

def query_pm_data(first_year, last_year, aqs, aqs_session, csa_code=None, cbsa_code=None, edt_ids=[0,5], par_codes= ['88101'],
                      pol_std_ids = [22], duration_codes=['V']):
    """Takes in a csa code, the first and last years, aqs, and session objects and returns a query of
    daily site record PM2.5 data for continuous and non-continuous methods

    Args:
        csa_code: an OMB combined statistical area code for the area to be analyzed
        cbsa_code: an OMB core based statistical area code.  Either a CSA code
                    or a CBSA code should be supplied, but not both.
        first_year: a four digit interger of the most recent year to analyze data
        last_year: a four digit interger of the most recent year to analyze data
        aqs: an SQLSoup object mapped to AQS
        aqs_session: an SQLAlchemy session object for AQS
        edt_ids: list of exceptional data type ids to use (default = 0 (no events)
                 and 5 (exclude regionally concurred events)
        par_codes: list of AQS parameter codes to pull data for
        pol_std_ids: list of ints for pollutant standard ids. default: 18 (2006 PM2.5 Annual Standard)
        duration_codes: list of strings of duration codes. (default: 'V' (Combined 24-hr))
    
    Returns: An SQLAlchemy Query object for PM2.5 data in the selected area

    Raises: None
    """
    first_date = datetime.date(first_year, 1, 1)
    last_date = datetime.date(last_year,12,31)

    query = aqs_session.query(aqs.site_daily_values.si_id)
    
    # Perform Joins
    query = query.join(aqs.site_basic, aqs.site_daily_values.si_id == aqs.site_basic.si_id)
    query = query.join(aqs.states, aqs.states.state_code == aqs.site_basic.state_code)

    # Join to CSA table or CBSA table
    query = query.join(aqs.counties, and_(aqs.site_basic.state_code == aqs.counties.stt_state_code,
                                                        aqs.site_basic.county_code == aqs.counties.county_code))
    query = query.outerjoin(aqs.cbsa_counties).outerjoin(aqs.core_based_statistical_areas)
    query = query.outerjoin(aqs.csa_cbsas).outerjoin(aqs.combined_statistical_areas)
    
        

    # Apply filters
    query = query.filter(and_(aqs.site_daily_values.sample_day >= first_date,
                              aqs.site_daily_values.sample_day <= last_date))
    query = query.filter(aqs.site_daily_values.edt_id.in_(edt_ids))
    query = query.filter(aqs.site_daily_values.parameter_code.in_(par_codes))
    query = query.filter(aqs.site_daily_values.pollutant_standard_id.in_(pol_std_ids))
    query = query.filter(aqs.site_daily_values.duration_code.in_(duration_codes))
    if csa_code != None:
        query = query.filter(aqs.combined_statistical_areas.csa_code == csa_code)
    elif cbsa_code != None:
        query = query.filter(aqs.core_based_statistical_areas.cbsa_code == cbsa_code)
    else:
        print 'No geographic filter applied.'

    site_code = (aqs.site_basic.state_code + '-' +
                    aqs.site_basic.county_code + '-' +
                    aqs.site_basic.site_id)

    query = query.add_columns(site_code.label('site_code'),
                              aqs.site_daily_values.parameter_code,
                              aqs.site_daily_values.sample_day,
                              aqs.site_daily_values.daily_value,
                              aqs.site_daily_values.duration_code,
                              aqs.site_daily_values.scheduled_date,
                              aqs.site_daily_values.creditable_day_indicator)
    
    query = query.order_by(site_code)

    #query = query.limit(100) # testing code
    
    return query

def query_sample_calendar(first_year, last_year, aqs, aqs_session):

    first_date = datetime.date(first_year, 1, 1)
    last_date = datetime.date(last_year,12,31)

    query = aqs_session.query(aqs.aqs_dates.aqs_date)
    query = query.filter(and_(aqs.aqs_dates.aqs_date >= first_date,
                              aqs.aqs_dates.aqs_date <= last_date))
    query = query.add_columns(aqs.aqs_dates.nams_every_3rd_day_ind,
                              aqs.aqs_dates.nams_every_6th_day_ind)
    return query

def query_site_collection_freq(first_year, last_year, aqs, aqs_session, csa_code=None, cbsa_code=None, par_codes= ['88101']):

    first_date = datetime.date(first_year, 1, 1)
    last_date = datetime.date(last_year,12,31)

    query = aqs_session.query(aqs.primary_monitor_periods.mo_id)
    
    # Perform Joins
    query = query.join(aqs.monitors, aqs.monitors.mo_id == aqs.primary_monitor_periods.mo_id)
    query = query.join(aqs.site_basic, aqs.site_basic.si_id == aqs.monitors.si_si_id)
    query = query.join(aqs.states, aqs.states.state_code == aqs.site_basic.state_code)
    query = query.join(aqs.req_coll_frequencies, aqs.primary_monitor_periods.mo_id == aqs.req_coll_frequencies.mo_mo_id) # Doesn't quite work, needs some filter on dates
    
    # Join to CSA table
    query = query.join(aqs.counties, and_(aqs.site_basic.state_code == aqs.counties.stt_state_code,
                                                        aqs.site_basic.county_code == aqs.counties.county_code))
    query = query.outerjoin(aqs.cbsa_counties).outerjoin(aqs.core_based_statistical_areas)
    query = query.outerjoin(aqs.csa_cbsas).outerjoin(aqs.combined_statistical_areas)

    # Apply filters
    query = query.filter(or_(aqs.primary_monitor_periods.end_date >= first_date,
                             aqs.primary_monitor_periods.end_date==None))
    query = query.filter(or_(aqs.req_coll_frequencies.req_coll_freq_end_date >= first_date,
                             aqs.req_coll_frequencies.req_coll_freq_end_date==None))
    query = query.filter(aqs.monitors.pa_parameter_code.in_(par_codes))

    if csa_code != None:
        query = query.filter(aqs.combined_statistical_areas.csa_code == csa_code)
    elif cbsa_code != None:
        query = query.filter(aqs.core_based_statistical_areas.cbsa_code == cbsa_code)
    else:
        print 'No geographic filter applied.'

    # Add columns
    site_code = (aqs.site_basic.state_code + '-' +
                    aqs.site_basic.county_code + '-' +
                    aqs.site_basic.site_id)

    monitor_code = (aqs.site_basic.state_code + '-' +
                    aqs.site_basic.county_code + '-' +
                    aqs.site_basic.site_id + '-' +
                    aqs.monitors.pa_parameter_code + '-' +
                    cast(aqs.monitors.poc, String(2)))

    query = query.add_columns(site_code.label('site_code'),
                              monitor_code.label('primary_monitor'),
                              aqs.req_coll_frequencies.cf_coll_freq_code,
                              aqs.primary_monitor_periods.begin_date.label('pri_mon_begin_date'),
                              aqs.primary_monitor_periods.end_date.label('pri_mon_end_date'),
                              aqs.req_coll_frequencies.req_coll_freq_begin_date,
                              aqs.req_coll_frequencies.req_coll_freq_end_date)
    return query
                              

print '-'*10, 'PM2.5 Imputation and Bootstrapping Tool', '-'*10
csa_or_cbsa = input("Please select the type of area to analyze:\n\
(1) A Combined Statistical Area (CSA)\n\
(2) A Core Based Statistical Area (CBSA)\n")
if csa_or_cbsa == 1:
    csa_code = raw_input('Please input the CSA code for the area you would like to analyze.\n')
    cbsa_code = None
elif csa_or_cbsa == 2:
    cbsa_code = raw_input('Please input the CBSA code for the area you would like to analyze.\n')
    csa_code = None
else: pass
first_year = input('Please input first year of the period you would like to analyze. \
(Typically use the most recent 5 year period)\n')
last_year = input('Please input last year of the period you would like to analyze.\n')

# Create New Directory
today = str(datetime.date.today())
orig_dir = os.getcwd()
if csa_code != None:
    geo_string = 'csa %s' %(csa_code)
elif cbsa_code != None:
    geo_string = 'cbsa %s' %(cbsa_code)
else:
    print 'No cbsa or csa code provided.'
directory = 'Bootstrapping Results %s %s' %(geo_string, today) 
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

raw_fname = 'PM25 Daily Site Values_%s_%s.txt' %(geo_string, today)
calendar_fname = 'aqs_sample_calendar.txt'
site_coll_freq_fname = 'site_collection_frequencies_%s.txt' %(today)


# Connect to AQS
aqs = mtss.connect_aqs()
aqs_session = aqs.session()

pm_query = query_pm_data(first_year, last_year, aqs, aqs_session, csa_code, cbsa_code)
sample_calendar_query = query_sample_calendar(first_year, last_year, aqs, aqs_session)
site_collection_freq_query = query_site_collection_freq(first_year, last_year, aqs, aqs_session, csa_code, cbsa_code)

print 'AQS Queries Built...'

mtss.write_query(pm_query, raw_fname)
mtss.write_query(sample_calendar_query, calendar_fname)
mtss.write_query(site_collection_freq_query, site_coll_freq_fname)

print "Done!"
endtime = datetime.datetime.now()
runtime = endtime - starttime
print "Run Time: %s" %(runtime)






