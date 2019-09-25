################################################################################
# Name : read_file_test
# Purpose : code to imitate how readcol in idl works for the mosaiq script
#Author : Benjamin Vaughan
#Start Date : September 25, 2019
#Additional Info
#
################################################################################
import numpy as np

def test(filename):
    inum, ramin, ramax, raavg, decmin, decmax, decavg, medianval, noise_key = np.loadtxt(filename, unpack=True)
    print(inum)
    print(ramin)
    print(ramax)
    print(raavg)
    print(decmin)
    print(decmax)
    print(decavg)
    print(medianval)
    print(noise_key)


if __name__ == '__main__':
    test('info_issa_map4.txt')
