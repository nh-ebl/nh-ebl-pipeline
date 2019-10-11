################################################################################
# Name : config.py
# Purpose : The purpose of this is to create a file that makes it easier to
# set up this code on a new computer
#Author : Benjamin Vaughan
#Start Date : Oct 11, 2019
#Additional Info
#
################################################################################

IrisLookupFile = 'info_issa_map4.txt' #please point to this textfile
# this means setting the path if it is in a different direct than the python
# scripts
IrisDir =  '../../../../IRISNOHOLES_B4H0' #please point to this
# to the directory where you are keeping all of your IRIS fits files.
fields = np.array([[177.4586, 32.6273],
                   [196.0334, 23.9452],
                   [161.8978, -26.7821],
                   [177.5121, 32.6308],
                   [196.0075, 23.9506],
                   [346.1128, -7.166],
                   [1.8064, -1.2497],
                   [270.6523, -14.6237],
                   [237.5430, -33.8045],
                   [268.4866, -34.7916],
                   [269.6507, -15.4847],
                   [231.3464, -17.9611],
                   [269.7189, -15.4893],
                   [237.3581, -33.7934]])
# these are the potential fields to look at, at some point it may be worthwile
# to have these written in a txtfile or something of that nature and then just
# read them into the config file with a call to a function such as np.loadtxt()
# for now, they are still hardcoded into the script
# index = 0 corresponds to the first field, index = 1 to the second, etc.
