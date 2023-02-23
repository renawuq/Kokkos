#ifndef UTILS_H_
#define UTILS_H_


template <class T>   
T FindMax(T a[],int n) //function to find largest element
{
	int i;
	T max;
	max = a[0];//assume that first element is max
	for(i=1;i<n;i++)
	{
		if(a[i]>max) //if currentelement is greater than max
			max =a[i]; //assign that number as max now

	}
	return max; //returns the largest number to main function
}

template <class T>  
T FindMin(T a[],int n) //function to find smallest element
{	
	int i;
	T min;
	min = a[0];// assuming first element as minimum
	for(i=1;i<n;i++)
	{
		if(a[i]<min)// If current element is smaller than min
			min =a[i];//assigning the smaller number to min
	}
	return min; //returns the smallest number to main function
}


// make the compiler ignore an unused variable
#ifndef UNUSED
#define UNUSED(x) ((void)(x))
#endif





#endif //	UTILS_H_