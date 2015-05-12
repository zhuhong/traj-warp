#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;

int main()
{
	cout<<"TEST"<<endl;
	vector<int> myvector;
	for (int i = 0; i < 6; ++i)
	{
		myvector.push_back(i*2);
	}
	for (vector<int>::iterator it = myvector.begin(); it != myvector.end(); it++)	
	{
		cout <<"\t"<<*it;
	}
	cout<<endl;

	
	int tv=3;
	// cout << "input a value"<<endl;
	// cin>>tv;
	vector<int>::iterator it;
	it = find(myvector.begin(), myvector.end(), tv);
	if(it != myvector.end()){
		cout<<"Element found..."<<*it<<endl;
	}else{
		cout<<"Element not found..."<<endl;	
	}
	
	cout<<fmod(-1.9, 2.0)+2<<endl;


	for (int i = 0; i < 6; ++i)
	{
		cout<<"i=\t"<<i<<endl;
		for (int j = 0; j < 6; ++j)
		{
			cout<<"j=\t"<<j<<"\t";

			if(j==3)
			{
				break;
			}
		}
		cout<<endl;
	}

}