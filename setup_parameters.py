path2bin = './bin/' # Set path to your environment executables

webservice = "IRIS"
network = "ZA" # NoMelt
compstr = "BH?"

# Start and end time to search for station list
tstart = "2011-12-08T00:00:00"
tend = "2012-12-09T00:00:00"

minmagnitude = 7.0 #6.5; 7.0; # Lower magnitude cutoff for event search    
search_dir = './' + webservice + '_' + network + '_' + str(minmagnitude) + '/' # OUTPUT PATH

input_stalist = 0 # 0 if use all stations
if input_stalist: # List of stations
    stalist = '/Users/russell/Lamont/ENAM/DATA/stalist_good.txt'
    text_file = open(stalist, "r")
    stations = text_file.read().split('\n')
    text_file.close()
    stations = ','.join(stations).replace(" ", "")
else: # Use all available stations
    stations = "*"

########## 1-run_BNGorient ##########
mindist = 30. # minimum distance to consider

########## 2-run_DLorient ##########

