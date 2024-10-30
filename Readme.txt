The script.m file has the code for decoding the data signals using the Weiner method and the Kalman filter
In the first line you may choose between DataSet1.mat, DataSet2.mat and DataSet3.mat
Within each datasat, in the second line you can choose between data.NoNoise_RxSignal, data.LowNoise_RxSignal and data.HighNoise_RxSignal

In the section starting in line 124 you can choose what method you want to use to decode the signal. 
In order to select a method you must uncomment the two lines for that method and comment the 2 lines for the other 2 methods.
You may choose between doing no channel estimation, Kalman filtering and Wiener FIlter

The script plots 2 graphs and displays an image of the decoded signal.


