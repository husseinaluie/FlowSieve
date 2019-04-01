
user="bstorer"
password="enter password here. Remove afterwards!"

motu_path="/Users/bstorer/Library/Python/3.7/lib/python/site-packages/motuclient.py"
motu_address="http://nrt.cmems-du.eu/motu-web/Motu"

LON_MIN="-80"
LON_MAX="-60"

LAT_MIN="25"
LAT_MAX="45"

date_min="2018-01-01 12:00:00"
date_max="2018-04-01 12:00:00"

depth_min="0.493"
depth_max="0.4942"

out_dir=`pwd`
out_name="hourly_data.nc"

# Daily average data
service_id="GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS"
product_id="global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh"

# Get the velocity and SSH
python ${motu_path}\
    --quiet \
    --user ${user}\
    --pwd ${password}\
    --motu ${motu_address}\
    --service-id ${service_id}\
    --product-id ${product_id}\
    --longitude-min ${LON_MIN}\
    --longitude-max ${LON_MAX}\
    --latitude-min ${LAT_MIN}\
    --latitude-max ${LAT_MAX}\
    --date-min ${date_min}\
    --date-max ${date_max}\
    --depth-min ${depth_min}\
    --depth-max ${depth_max}\
    --variable uo\
    --variable vo\
    --variable zos\
    --variable thetao\
    --out-dir ${out_dir}\
    --out-name ${out_name}

