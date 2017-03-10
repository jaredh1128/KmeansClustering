#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

double** GetEuclideanDistanceMatrix(double**, double**, int, int, int);
double** GetManhattanDistanceMatrix(double**, double**, int, int, int);
double** GetMeans(vector<string>, double**, double**, int, int, int, int*);
double** GetData(vector<string>, int);
vector<string> GetDimensions();

void Euclidean(vector<string> titles, double** dataset, int iclusters, int datasize, int dimensionsize) //performs Euclidean distance clustering
{
	double** distancematrix = new double*[iclusters];
	double** clustercenters = new double*[iclusters];
	double** clustercompare = new double*[iclusters];
	int nochange = 1;
	int index = 0;
	int compare;
	int iteration = 0;
	int counter = 0;
	int* clusterindex = new int[datasize];

	for (int i = 0; i < iclusters; i++)
		distancematrix[i] = new double[datasize];

	for (int i = 0; i < iclusters; i++)
	{
		clustercenters[i] = new double[dimensionsize];
		clustercompare[i] = new double[dimensionsize];
	}
	cout << "Performing K-Means algorithm through Euclidean distance" << endl << "-------------------------------------------------------" << endl << endl;
	cout << "Original Cluster Centers: " << endl;
	for (int t = 0; t < dimensionsize; t++)
	{
		cout << setw(10) << titles.at(t);
	}
	cout << endl;
	for (int i = 0; i < iclusters; i++)
	{
		for (int j = 0; j < dimensionsize; j++)
		{
			clustercenters[i][j] = dataset[i][j]; //sets the first k clusters as centroids
			clustercompare[i][j] = clustercenters[i][j];
			cout << setw(10) << clustercenters[i][j];
		}
		cout << endl;
	}
	while (nochange != 0) //if a single point differs from the previous iteration, perform another iteration
	{
		nochange = 0;
		iteration++;
		cout << endl << "Iteration " << iteration << ": " << endl << "----------" << endl;
		distancematrix = GetEuclideanDistanceMatrix(dataset, clustercenters, iclusters, datasize, dimensionsize);

		cout << "Distance Matrix Results: " << endl;
		for (int i = 0; i < iclusters; i++)
		{
			cout << "Distance from Cluster " << i + 1 << ": ";
			for (int j = 0; j < datasize; j++)
			{
				cout << distancematrix[i][j] << ", ";
			}
			cout << endl;
		}
		cout << "Groups: ";
		for (int i = 0; i < datasize; i++)
		{
			for (int j = 0; j < iclusters - 1; j++)
			{
				if (counter == 0)
				{
					compare = distancematrix[j][i];
					index = 0;
					counter++;
				}
				if (compare > distancematrix[j + 1][i])
				{
					index = j + 1;
					compare = distancematrix[j + 1][i];
				}
			}
			counter = 0;
			clusterindex[i] = index; //Array of cluster indices (acts as G from the notes)
			cout << index << " ";
		}
		cout << endl;

		clustercompare = clustercenters; //Create a copy of current cluster centers to compare at a later time
		clustercenters = GetMeans(titles, dataset, clustercenters, iclusters, datasize, dimensionsize, clusterindex);
		cout << "New Cluster Centers: " << endl;
		for (int t = 0; t < dimensionsize; t++)
		{
			cout << setw(10) << titles.at(t);
		}
		cout << endl;
		for (int i = 0; i < iclusters; i++)
		{
			for (int j = 0; j < dimensionsize; j++)
			{
				if (clustercenters[i][j] != clustercompare[i][j]) //check if the current and previous iterations have changed
					nochange++;
				cout << setw(10) << clustercenters[i][j];
			}
			cout << endl;
		}
		cout << endl << "Old Cluster Centers: " << endl;
		for (int t = 0; t < dimensionsize; t++)
		{
			cout << setw(10) << titles.at(t);
		}
		cout << endl;
		for (int i = 0; i < iclusters; i++)
		{
			for (int j = 0; j < dimensionsize; j++)
			{
				cout << setw(10) << clustercompare[i][j];
			}
			cout << endl;
		}
	}
	cout << "Done" << endl << endl;
}

void Manhattan(vector<string> titles, double** dataset, int iclusters, int datasize, int dimensionsize) //performs Manhattan distance clustering
{
	double** distancematrix = new double*[iclusters];
	double** clustercenters = new double*[iclusters];
	double** clustercompare = new double*[iclusters];
	int nochange = 1;
	int index = 0;
	int compare;
	int iteration = 0;
	int counter = 0;
	int* clusterindex = new int[datasize];

	for (int i = 0; i < iclusters; i++)
		distancematrix[i] = new double[datasize];

	for (int i = 0; i < iclusters; i++)
	{
		clustercenters[i] = new double[dimensionsize];
		clustercompare[i] = new double[dimensionsize];
	}
	cout << "Performing K-Means algorithm through Manhattan distance" << endl << "-------------------------------------------------------" << endl << endl;
	cout << "Original Cluster Centers: " << endl;
	for (int t = 0; t < dimensionsize; t++)
	{
		cout << setw(10) << titles.at(t);
	}
	cout << endl;
	for (int i = 0; i < iclusters; i++)
	{
		for (int j = 0; j < dimensionsize; j++)
		{
			clustercenters[i][j] = dataset[i][j]; //sets the first k clusters as centroids
			clustercompare[i][j] = clustercenters[i][j];
			cout << setw(10) << clustercenters[i][j];
		}
		cout << endl;
	}
	while (nochange != 0) //if a single point differs from the previous iteration, perform another iteration
	{
		nochange = 0;
		iteration++;
		cout << endl << "Iteration " << iteration << ": " << endl << "----------" << endl;
		distancematrix = GetManhattanDistanceMatrix(dataset, clustercenters, iclusters, datasize, dimensionsize);

		cout << "Distance Matrix Results: " << endl;
		for (int i = 0; i < iclusters; i++)
		{
			cout << "Distance from Cluster " << i + 1 << ": ";
			for (int j = 0; j < datasize; j++)
			{
				cout << distancematrix[i][j] << ", ";
			}
			cout << endl;
		}
		cout << "Groups: ";
		for (int i = 0; i < datasize; i++)
		{
			for (int j = 0; j < iclusters - 1; j++)
			{
				if (counter == 0)
				{
					compare = distancematrix[j][i];
					index = 0;
					counter++;
				}
				if (compare > distancematrix[j + 1][i])
				{
					index = j + 1;
					compare = distancematrix[j + 1][i];
				}
			}
			counter = 0;
			clusterindex[i] = index; //Array of cluster indices (acts as G from the notes)
			cout << index << " ";
		}
		cout << endl;

		clustercompare = clustercenters; //Create a copy of current cluster centers to compare at a later time
		clustercenters = GetMeans(titles, dataset, clustercenters, iclusters, datasize, dimensionsize, clusterindex);
		cout << "New Cluster Centers: " << endl;
		for (int t = 0; t < dimensionsize; t++)
		{
			cout << setw(10) << titles.at(t);
		}
		cout << endl;
		for (int i = 0; i < iclusters; i++)
		{
			for (int j = 0; j < dimensionsize; j++)
			{
				if (clustercenters[i][j] != clustercompare[i][j]) //check if the current and previous iterations have changed
					nochange++;
				cout << setw(10) << clustercenters[i][j];
			}
			cout << endl;
		}
		cout << endl << "Old Cluster Centers: " << endl;
		for (int t = 0; t < dimensionsize; t++)
		{
			cout << setw(10) << titles.at(t);
		}
		cout << endl;
		for (int i = 0; i < iclusters; i++)
		{
			for (int j = 0; j < dimensionsize; j++)
			{
				cout << setw(10) << clustercompare[i][j];
			}
			cout << endl;
		}
	}
	cout << "Done" << endl << endl;
}
double** GetMeans(vector<string> titles, double** dataset, double** clustercenters, int iclusters, int datasize, int dimensionsize, int* clusterindex)
{
	double** newclustercenters = new double*[iclusters];
	double*** pointdisplay = new double**[iclusters];
	int* counter = new int[iclusters];

	for (int i = 0; i < iclusters; i++)
	{
		pointdisplay[i] = new double*[datasize];
		for (int j = 0; j < datasize; j++)
		{
			pointdisplay[i][j] = new double[dimensionsize];
		}
	}

	for (int i = 0; i < iclusters; i++)
	{
		newclustercenters[i] = new double[dimensionsize];
	}
	for (int i = 0; i < iclusters; i++)
	{
		counter[i] = 0;

		for (int j = 0; j < dimensionsize; j++)
		{
			newclustercenters[i][j] = 0;
		}
	}
	for (int j = 0; j < datasize; j++)
	{
		for (int k = 0; k < dimensionsize; k++)
		{
			newclustercenters[clusterindex[j]][k] += dataset[j][k];//Add datapoints under a specific cluster together
		}
		counter[clusterindex[j]]++; //increment the count of the cluster that the datapoint belongs to
	}

	for (int j = 0; j < datasize; j++)
	{
		for (int k = 0; k < dimensionsize; k++)
		{
			pointdisplay[clusterindex[j]][j][k] = dataset[j][k];
		}
	}

	for (int i = 0; i < iclusters; i++)
	{
		cout << "Cluster " << i + 1 << ": " << endl;
		for (int t = 0; t < dimensionsize; t++)
		{
			cout << setw(10) << titles.at(t);
		}
		cout << endl;
		for (int j = 0; j < datasize; j++)
		{
			for (int k = 0; k < dimensionsize; k++)
			{
				if (pointdisplay[i][j][k] > -6e66)
					cout << setw(10) << pointdisplay[i][j][k];
				if (pointdisplay[i][j][k] > -6e66 && k == dimensionsize - 1)
					cout << endl;
			}
		}
		cout << "-------------------" << endl;
	}
	for (int i = 0; i < iclusters; i++)
	{
		for (int j = 0; j < dimensionsize; j++)
		{
			newclustercenters[i][j] = newclustercenters[i][j] / counter[i]; //take the mean of the point sums by dividing by the cluster counters
		}
	}
	return newclustercenters;
}

double** GetEuclideanDistanceMatrix(double** dataset, double** clustercenters, int iclusters, int datasize, int dimensionsize)
{
	double** distances = new double*[iclusters];
	double calcdist = 0;
	for (int i = 0; i < iclusters; i++)
	{
		distances[i] = new double[datasize];
	}
	for (int k = 0; k < iclusters; k++)
	{
		for (int i = 0; i < datasize; i++)
		{
			for (int j = 0; j < dimensionsize; j++)
			{
				calcdist += pow((dataset[i][j] - clustercenters[k][j]), 2); //continue adding squares of point differences for n-dimensions
			}
			calcdist = sqrt(calcdist);
			distances[k][i] = calcdist;
			calcdist = 0;
		}
	}
	return distances;
}

double** GetManhattanDistanceMatrix(double** dataset, double** clustercenters, int iclusters, int datasize, int dimensionsize)
{
	double** distances = new double*[iclusters];
	double calcdist = 0;
	for (int i = 0; i < iclusters; i++)
	{
		distances[i] = new double[datasize];
	}
	for (int k = 0; k < iclusters; k++)
	{
		for (int i = 0; i < datasize; i++)
		{
			for (int j = 0; j < dimensionsize; j++)
			{
				calcdist += abs(dataset[i][j] - clustercenters[k][j]); //continue adding absolute values of point differences for n-dimensions
			}
			distances[k][i] = calcdist;
			calcdist = 0;
		}
	}
	return distances;
}
double** GetData(vector<string> titles, int idatasize)
{
	double** datasets;
	double data;

	datasets = new double*[idatasize];
	for (int i = 0; i < idatasize; i++)
	{
		datasets[i] = new double[titles.size()];
	}

	for (int i = 0; i < idatasize; i++)
	{
		for (int j = 0; j < titles.size(); j++)
		{
			cout << "Please enter " << titles.at(j) << " data point " << i << ":" << endl;
			cin >> data;
			datasets[i][j] = data; //enter data from left-right/top-bottom formation
			while (cin.fail()) //input validation
			{
				cout << "Invalid entry" << endl;
				cout << "Please enter " << titles.at(i) << " data point " << j << ":" << endl;
				cin.clear();
				cin.ignore(256, '\n');
				cin >> data;
				datasets[i][j] = data;
			}
		}
	}
	return datasets;
}
vector<string> GetDimensions()
{
	int dimensions;
	vector<string> titlestorage;
	string title;

	cout << "Enter number of dimensions: " << endl;
	cin >> dimensions;
	while (cin.fail())//input validation
	{
		cout << "Invalid input " << endl;
		cout << "Enter number of dimensions: " << endl;
		cin.clear();
		cin.ignore(256, '\n');
		cin >> dimensions;
	}

	for (int i = 0; i < dimensions; i++)
	{
		cout << "Please enter title " << i + 1 << ": " << endl;
		cin >> title;
		titlestorage.push_back(title);
	}
	return titlestorage;
}

void main()
{
	double** datapoints;
	int clusters;
	vector<string> titlestorage;
	int dimensionsize;
	int datasize;

	cout << "Enter number of desired clusters: " << endl;
	cin >> clusters;
	while (cin.fail() || clusters <= 0)
	{
		cout << "Invalid entry" << endl;
		cout << "Enter number of desired clusters: " << endl;
		cin.clear(); //clear input
		cin.ignore(256, '\n');	//ignore error bit
		cin >> clusters;
	}

	titlestorage = GetDimensions(); //Gets titlenames and dimension size
	cout << "Enter maximum number of data points per dimension: " << endl;
	cin >> datasize;
	while (cin.fail() || datasize < clusters)
	{
		cout << "Invalid entry" << endl;
		cout << "Enter maximum number of data points per dimension: " << endl;
		cin.clear();
		cin.ignore(256, '\n');
		cin >> datasize;
	}
	datapoints = new double*[datasize];
	for (int i = 0; i < datasize; i++)
	{
		datapoints[i] = new double[titlestorage.size()];
	}
	datapoints = GetData(titlestorage, datasize); //Creates a datasize X dimensionsize matrix containing all input data

	Euclidean(titlestorage, datapoints, clusters, datasize, titlestorage.size()); //performs clustering using Euclidean distance calculation
	Manhattan(titlestorage, datapoints, clusters, datasize, titlestorage.size()); //performs clustering using Manhattan distance calculation

	cout << "Demo Done" << endl;
	system("pause");
}