# script to convert a CSV from URA or Ellipsometry measurement to a 2 column csv
# use line count if statement to ignore a number of initial rows
# use row argument to select rows to convert
# ensure yData is from 0 to 100 not 0 to 1
# Author: Sol Cotton

from __future__ import print_function
import csv

F = open("Edmund_450_60.csv", "w+")
with open('rp_f_7_60deg.csv') as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=',')
	line_count = 0
	for row in csv_reader:
		if line_count==0:
			line_count+=1
		else:
			writeline = str(row[0]) + "," + str(row[1]) + "\n"
			F.write(writeline)
			line_count += 1
F.close()