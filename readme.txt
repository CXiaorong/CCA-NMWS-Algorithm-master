Operating environment：Windows 10，R3.4.1

Input datas:
	(1)binary mutation matrix.   data\GBM\SNVdata_440.csv     data\OVCA\SNVdata_2547.csv   data\THCA\SNVdata_3420.csv
	(2)network matrix.	     data\GBM\network_440.csv    data\OVCA\network_2547.csv   data\THCA\SNVdata_3420.csv
			
Steps:
	1.First select the function(1)-function(14) function and run it.
	         function(1): Calculate the actual number of edges in the network
	         function(2): fitness function--NMWS model
	         function(3): Compute the selection probability function.
	         function(4): crossover function 
	         function(5): Mutation function
	         function(6): mutation eg: gene n=1,2,3 ->n=1,2,6
	         function(7): select function
	         function(8): competition function
	         function(9): Compare the magnitude of two numbers
	         function(10): GA2 function
	         function(11): GA4 function
	         function(12): Cooperative pool evolution process.
	         function(13): Initialize the population function
	         function(14): significance test function
	2.load data
	3.Initialization parameters.
	4.Iterative loop

Output:
	A gene set and its corresponding fitness function value.
	For example: GBM dataset, k=5, the result is: "CDKN2A" "MDM2"   "MDM4"   "RB1"    "TP53"   "1.05083333333333"
