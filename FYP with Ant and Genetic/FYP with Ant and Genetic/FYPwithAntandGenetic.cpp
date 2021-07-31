// FYPwithAntandGenetic.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <functional> 
#include <iterator>
#include <math.h>
#include <list>
#include <chrono>
using std::cout;
using std::endl;
using std::array;
using std::vector;


#define numberOfCities 10
#define maxDistanceBetweenCities 50
#define maxRouteLength (maxDistanceBetweenCities*numberOfCities)
#define numberOfAnts 10
#define numberOfGenes 10

//algorithm constants
//heuristic functions 
#define ALPHA 1.0 //importance of trail
#define BETA 5.0 //importance of visibility
#define evaporationRateFactor 0.5 
#define QVal 100
#define maxPaths 20
#define numberOfGenerations 20
#define maxCombinations (maxPaths * numberOfCities)
#define initial_Pheromones (1.0/numberOfCities) //weight of the edge


struct cityLocation
{
	int x, y;
};

struct GenesType
{
	int gene[numberOfCities];
	int fitness;
};

struct AntsType
{
	int currentCity, nextCity, pathNumhber;
	int citiesVisted[numberOfCities];
	int pathTaken[numberOfCities];
	float pathDistance;
};


cityLocation cities[numberOfCities];
AntsType ants[numberOfAnts];
GenesType bestOrder;
vector<GenesType> startingGenes;
vector<GenesType> genes;
vector<GenesType> newGenes;
vector<GenesType> eliteGenes;
vector<GenesType> crossOverGenes;
vector<GenesType> newRandomGenes;
float fitness[numberOfGenes];
double distance[numberOfCities][numberOfCities];
double pheromones[numberOfCities][numberOfCities];
double best = (double)maxRouteLength;
double bestDistance = (double)maxRouteLength;
double record = (double)maxRouteLength;
int bestPathNumber = NULL;
int bestPath[6];
bool Hybrid = false;
double tourTotal;
double tourAverage;
double cycleTotal;
double cycleAverage;



void calculateRouteDistance()
{
	int cityFrom, cityTo;
	cycleTotal= 0;
	cycleAverage= 0;

	for (int i = 0; i < numberOfGenes; i++)
	{
		GenesType currentOrder;
		currentOrder = genes[i];

		//float disToNextCity = 0;
		float routeDistance = 0;
		//cout << endl;
		//cout << endl;
		for (int j = 0; j < numberOfCities-1; j++)
		{
			//cout << j << endl;
			cityFrom = currentOrder.gene[j];
			cityTo = currentOrder.gene[j+1];
			/*cout << "city From" << cityFrom << endl;
			//cout << "city to" << cityTo << endl;
			//cout << "distance from (" << cityFrom << ") to (" << cityTo << ") :" << distance[cityFrom][cityTo] << endl;
			*/
			routeDistance += distance[cityFrom][cityTo];
			//cout << routeDistance << endl;
		}
		cityFrom = currentOrder.gene[numberOfCities - 1];
		cityTo = currentOrder.gene[0];
		//cout << "distance from (" << cityFrom << ") to (" << cityTo << ") :" << distance[cityFrom][cityTo] << endl;
		//cout << "distance from (" << cityFrom << ") to (" << cityTo << ") :" << distance[numberOfCities - 1][0] << endl;

		routeDistance += distance[cityFrom][cityTo];
		//cout << endl;
		/*for (int j = 0; j < numberOfCities; j++)
		{
			cout << genes[i].gene[j] << ", ";
		}
		cout << endl;
		cout << routeDistance << endl;
		*/
		if(i < numberOfGenes*0.8)
		{
			cycleTotal += routeDistance;
		}
		if (routeDistance < bestDistance)
		{
			bestDistance = routeDistance;
			bestOrder = genes[i];
		}
		fitness[i] = 1 / (routeDistance + 1);
		//cout << fitness[i] << endl;
	}
	cycleAverage = cycleTotal / numberOfGenes*0.8;
	cout << "cycle average: " << cycleAverage << endl;
	//cout << endl;
	float total = 0;
	for (int i = 0; i < sizeof(fitness) / sizeof(int); i++)
	{
		total += fitness[i];
		//cout << "size of fitness" << sizeof(fitness) / sizeof(int) << endl;
		//cout << fitness[i] << endl;
		//cout << total << endl;
	}
	for (int i = 0; i < sizeof(fitness) / sizeof(int); i++)
	{
		//cout << sizeof(fitness) / sizeof(int) << endl;
		fitness[i] = fitness[i] / total;
		//cout << fitness[i] << endl;
	}

	//cout << endl;
	//cout << endl;
	cout << "This is the best distance: " << bestDistance << endl;
	/*cout << "This is the best route: ";
	for (int i = 0; i < numberOfCities; i++)
	{
		cout << bestOrder.gene[i] << ", ";
	}
	cout << endl;*/

}


void generateGenes()
{
	//startingGenes.clear();
	int initialorder[numberOfCities];
	int order[numberOfCities];
	//create initial order
	for (int i = 0; i < numberOfCities; i++)
	{
		initialorder[i] = i;
	}


	for (int i = 0; i < numberOfGenes; i++)
	{
		GenesType NewGene;
		for (int j = 0; j < numberOfCities; j++)
		{

			NewGene.gene[j] = initialorder[j];
		}

		std::random_shuffle(std::begin(NewGene.gene), std::end(NewGene.gene));
		startingGenes.push_back(NewGene);
/*
		for (int j = 0; j < numberOfCities; j++)
		{
			cout << startingGenes[i].gene[j] << ",";
		}
		cout << endl;
		*/
	}
}

double antInformation(int cityFrom, int cityTo);
GenesType mutate(GenesType & geneToMutate)
{
	int cityFrom, cityTo;
	int worstTo = NULL;
	int secondWorstTo= NULL;
	int worstToPosition = NULL;
	int secondWorstToPosition = NULL;
	float worstSegment = INFINITY;
	float secondWorstSegment = INFINITY;
	double equationDenominator = 0.0;
	/*/
	for (int i = 0; i < numberOfCities; i++)
	{
		cout << geneToMutate.gene[i] << ", ";
	}
	cout << endl;
	*/
	
	for (int i = 0; i < numberOfCities - 1; i++)
	{
		cityFrom = geneToMutate.gene[i];
		cityTo = geneToMutate.gene[i + 1];
		antInformation(cityFrom, cityTo);
		
		equationDenominator += antInformation(cityFrom, cityTo);
		//cout << antInformation(cityFrom, cityTo) << endl;
		//cout << equationDenominator << endl;
	}
	for (int i = 0; i < numberOfCities - 1; i++)
	{
		cityFrom = geneToMutate.gene[i];
		cityTo = geneToMutate.gene[i + 1];
		double information = antInformation(cityFrom, cityTo) / equationDenominator;
		//cout << information << endl;
		if (information < worstSegment)
		{
			secondWorstSegment = worstSegment;
			worstSegment = information;
			secondWorstTo = worstTo;
			worstTo = cityTo;
			secondWorstToPosition = worstToPosition;
			worstToPosition = i+1;
			//cout <<"new worst"<< worstTo << endl;
			//cout << "new second worst" << secondWorstTo << endl;
		} else if (information < secondWorstSegment)
		{
			secondWorstSegment = information;
			secondWorstTo = cityTo;

			secondWorstToPosition = i + 1;
			//cout << "new second worst" << secondWorstTo << endl;
		}
	}

	//cout << "worst" << worstTo << endl;
	//cout << "second worst" << secondWorstTo << endl;
	std::swap(geneToMutate.gene[worstToPosition], geneToMutate.gene[secondWorstToPosition]);

	return geneToMutate;
}



GenesType CrossOverGeneraterWithMutation(vector<GenesType> & crossOverGenes)
{
	int *exists;
	int firstGene = rand() % crossOverGenes.size();
	int secondGene = rand() % crossOverGenes.size();
	/*	if (secondGene == firstGene)
	{
	secondGene += 1;
	}*/
	int newCrossOverStart = rand() % numberOfCities;
	GenesType newCrossOverGene;
	//cout << endl;
	//first part of need crossed over gene
	for (int i = 0; i < newCrossOverStart + 1; i++)
	{
		newCrossOverGene.gene[i] = crossOverGenes[firstGene].gene[i];
		//cout << newCrossOverGene.gene[i] << ",";
	}
	//remaining part of new gene
	for (int i = newCrossOverStart + 1; i < numberOfCities; i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
			int search = crossOverGenes[secondGene].gene[j];
			bool addValue = true;
			for (int k = 0; k < newCrossOverStart + 1; k++)
			{
				if (crossOverGenes[secondGene].gene[j] == crossOverGenes[firstGene].gene[k])
				{
					addValue = false;
				}
			}
			if (addValue == true)
			{
				newCrossOverGene.gene[i] = crossOverGenes[secondGene].gene[j];
				//	cout << newCrossOverGene.gene[i] << ",";
				i++;
			}

		}
	}

	GenesType mutatedGene = mutate(newCrossOverGene);
	newCrossOverGene = mutatedGene;
	/*
	for (int i = 0; i < numberOfCities; i++)
	{
	cout << newCrossOverGene.gene[i] << ", ";
	}
	cout << endl;

	//cout << endl;
	for (int i = 0; i < numberOfCities; i++)
	{
	//	cout << crossOverGenes[firstGene].gene[i] << ",";
	}
	//cout << endl;
	for (int i = 0; i < numberOfCities; i++)
	{
	//cout << crossOverGenes[secondGene].gene[i] << ",";
	}
	//cout << endl;
	for (int i = 0; i < numberOfCities; i++)
	{
	//cout << newCrossOverGene.gene[i] << ",";
	}
	//	cout << endl;
	*/


	return newCrossOverGene;
}



GenesType CrossOverGenerater(vector<GenesType> & crossOverGenes)
{
	int *exists;
	int firstGene = rand() % crossOverGenes.size();
	int secondGene = rand() % crossOverGenes.size();
	/*	if (secondGene == firstGene)
	{
	secondGene += 1;
	}*/
	int newCrossOverStart = rand() % numberOfCities;
	GenesType newCrossOverGene;
	//cout << endl;
	//first part of need crossed over gene
	for (int i = 0; i < newCrossOverStart + 1; i++)
	{
		newCrossOverGene.gene[i] = crossOverGenes[firstGene].gene[i];
		//cout << newCrossOverGene.gene[i] << ",";
	}
	//remaining part of new gene
	for (int i = newCrossOverStart + 1; i < numberOfCities; i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
			int search = crossOverGenes[secondGene].gene[j];
			bool addValue = true;
			for (int k = 0; k < newCrossOverStart + 1; k++)
			{
				if (crossOverGenes[secondGene].gene[j] == crossOverGenes[firstGene].gene[k])
				{
					addValue = false;
				}
			}
			if (addValue == true)
			{
				newCrossOverGene.gene[i] = crossOverGenes[secondGene].gene[j];
			//	cout << newCrossOverGene.gene[i] << ",";
				i++;
			}

		}
	}

	//GenesType mutatedGene = mutate(newCrossOverGene);
	//newCrossOverGene = mutatedGene;
	/*
	for (int i = 0; i < numberOfCities; i++)
	{
		cout << newCrossOverGene.gene[i] << ", ";
	}
	cout << endl;
	
	//cout << endl;
	for (int i = 0; i < numberOfCities; i++)
	{
	//	cout << crossOverGenes[firstGene].gene[i] << ",";
	}
	//cout << endl;
	for (int i = 0; i < numberOfCities; i++)
	{
		//cout << crossOverGenes[secondGene].gene[i] << ",";
	}
	//cout << endl;
	for (int i = 0; i < numberOfCities; i++)
	{
		//cout << newCrossOverGene.gene[i] << ",";
	}
//	cout << endl;
	*/


	return newCrossOverGene;
}



void selectGenes(vector<GenesType> & geneArray, float prob[numberOfGenes])
{

	float temp;
	array<int, 6> tempArray;
	vector<GenesType> newCrossOverArray;
	eliteGenes.clear();
	crossOverGenes.clear();
	newRandomGenes.clear();

	bool swapped = true;
	while (swapped == true)
	{
		swapped = false;
		for (int i = 0; i < (sizeof(fitness) / sizeof(int)) - 1; i++)
		{
			if (fitness[i] < fitness[i + 1])
			{
				std::swap(fitness[i], fitness[i + 1]);
				std::swap(geneArray[i], geneArray[i + 1]);
				swapped = true;
			}
		}
	}
	for (int i = 0; i < sizeof(fitness) / sizeof(int); i++)
	{
		//cout << fitness[i] << endl;
		if (i < numberOfGenes * 0.2)
		{
			eliteGenes.push_back(geneArray[i]);
		}
		if (i < numberOfGenes * 0.6)
		{
			crossOverGenes.push_back(geneArray[i]);
		}
	}
	if (eliteGenes.size() + crossOverGenes.size() < numberOfGenes)
	{
		int currentGenes = eliteGenes.size() + crossOverGenes.size();
		int randomNeeded = numberOfGenes - currentGenes;
		int initialorder[numberOfCities];
		for (int i = 0; i < numberOfCities; i++)
		{
			initialorder[i] = i;
		}
		for (int i = 0; i < randomNeeded; i++)
		{
			GenesType newRandomGene;
			for (int j = 0; j < numberOfCities; j++)
			{

				newRandomGene.gene[j] = initialorder[j];
			}
			std::random_shuffle(std::begin(newRandomGene.gene), std::end(newRandomGene.gene));
			newRandomGenes.push_back(newRandomGene);
		}
	}


	for (int i = 0; i < eliteGenes.size(); i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
			//cout << eliteGenes[i].gene[j] << ",";
		}
		//cout << endl;
	}
	//cout << endl;
	for (int i = 0; i < crossOverGenes.size(); i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
			//cout << crossOverGenes[i].gene[j] << ",";
		}
	//	cout << endl;
	}
//	cout << endl;
	for (int i = 0; i < newRandomGenes.size(); i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
			//cout << newRandomGenes[i].gene[j] << ",";
		}
	//	cout << endl;
	}
	//cout << endl;
	for (int i = 0; i < crossOverGenes.size(); i++)
	{
		GenesType newCrossedOverGene;
		if (Hybrid == false) 
		{
			newCrossedOverGene = CrossOverGenerater(crossOverGenes);
		}
		else
		{
			newCrossedOverGene = CrossOverGeneraterWithMutation(crossOverGenes);
		}
		
		newCrossOverArray.push_back(newCrossedOverGene);
	}
	//cout << crossOverGenes.size() << endl;;
	//cout << newCrossOverArray.size() << endl;
	crossOverGenes.clear();
	for (int i = 0; i < newCrossOverArray.size(); i++)
	{
		crossOverGenes.push_back(newCrossOverArray[i]);
	}
	newGenes.clear();


	for (int i = 0; i < eliteGenes.size(); i++)
	{
		newGenes.push_back(eliteGenes[i]);
		//cout << endl;
	}
	//cout << endl;
	for (int i = 0; i < crossOverGenes.size(); i++)
	{
		newGenes.push_back(crossOverGenes[i]);
		//cout << endl;
	}
	for (int i = 0; i < newRandomGenes.size(); i++)
	{
		newGenes.push_back(newRandomGenes[i]);
		//cout << endl;
	}

	for (int i = 0; i < numberOfGenes; i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
		//	cout << newGenes[i].gene[j] << ",";
		}
		//cout << endl;
	}
	genes.clear();
	for (int i = 0; i < numberOfGenes; i++)
	{
		genes.push_back(newGenes[i]);
	}

	//	return newGenes;

}



void initiate()
{
	int cityFrom, cityTo, ant;
	//cout << endl;
	//generate the cities
	for (cityFrom = 0; cityFrom < numberOfCities; cityFrom++)
	{
		cities[cityFrom].x = rand() % maxDistanceBetweenCities;
		cities[cityFrom].y = rand() % maxDistanceBetweenCities;
		cout << cities[cityFrom].x << "," << cities[cityFrom].y << endl;

		for (cityTo = 0; cityTo < numberOfCities; cityTo++)
		{
			distance[cityFrom][cityTo] = 0.0;
			pheromones[cityFrom][cityTo] = initial_Pheromones;
		}
	}

	//cout << endl;
	//cout << endl;
	//calculate the distance
	for (cityFrom = 0; cityFrom < numberOfCities; cityFrom++)
	{
		//cout << "city from " << cityFrom << "	";
		for (cityTo = 0; cityTo < numberOfCities; cityTo++)
		{
			//cout << "city to "<<  cityTo << " : ";
			if (cityTo != cityFrom && distance[cityFrom][cityTo] == 0.0)
			{
				int xDistance = pow(abs(cities[cityFrom].x - cities[cityTo].x), 2);
				int yDistance = pow(abs(cities[cityFrom].y - cities[cityTo].y), 2);

				distance[cityFrom][cityTo] = sqrt(xDistance + yDistance);
				distance[cityTo][cityFrom] = distance[cityFrom][cityTo];
				//cout << "distance from (" << cityFrom<< ") to (" << cityTo<< ") :" << distance[cityFrom][cityTo] << endl;
			}
			/*else
			{
			cout << ";	";
			}*/
		}
		//cout << endl;

	}


	cityTo = 0;
	cout << endl;
	cout << endl;
	for (ant = 0; ant < numberOfAnts; ant++)
	{
		//cout << "ant number: " << ant << endl;
		if (cityTo == numberOfCities)
		{
			cityTo = 0;
		}

		ants[ant].currentCity = cityTo++;
		for (cityFrom = 0; cityFrom < numberOfCities; cityFrom++)
		{
			ants[ant].citiesVisted[cityFrom] = 0;
			ants[ant].pathTaken[cityFrom] = -1;
		}

		ants[ant].pathNumhber = 1;
		ants[ant].pathTaken[0] = ants[ant].currentCity;
		ants[ant].nextCity = -1;
		ants[ant].pathDistance = 0;

		ants[ant].citiesVisted[ants[ant].currentCity] = 1;
	}


}


void restartAllAnts()
{
	int ant, cityTo = 0;
	tourTotal = 0;
	tourAverage = 0;

	for (ant = 0; ant < numberOfAnts; ant++)
	{
		//cout << best << endl;
		if(ant < numberOfAnts*0.8)
		{
			tourTotal += ants[ant].pathDistance;
		}
		//cout << ants[ant].pathDistance << endl;
		//cout << endl;
		//cout << "best distance: " << best << endl;
		//cout << "path distance: " << ants[ant].pathDistance << endl;
		if (ants[ant].pathDistance < best)
		{
			best = ants[ant].pathDistance;
			//cout << "best distance" << ants[ant].pathDistance << endl;
			bestPathNumber = ant;
			cout << "New best distance :" << best << endl;
		}
/*
		for (int i = 0; i < numberOfCities; i++)
		{
			cout << ants[ant].pathTaken[i] << ", ";
		}
		cout << endl;*/
		if (Hybrid == true)
		{
			
			AntsType currentAnt;
			GenesType NewGene;
			currentAnt = ants[ant];
		/*	cout << "current ant: ";
			for (int i = 0; i < numberOfCities; i++)
			{
				cout << ants[ant].pathTaken[i] << ", ";
			}
			cout << endl;*/
			//cout << "new gene: ";
			for (int i = 0; i < numberOfCities; i++)
			{
				NewGene.gene[i] = currentAnt.pathTaken[i];
				//cout << ants[ant].pathTaken[i] << ", ";
			}
		//	cout << endl;
			genes.push_back(NewGene);


		}


		ants[ant].nextCity = -1;
		ants[ant].pathDistance = 0.0;

		for (int i = 0; i < numberOfCities; i++)
		{
			ants[ant].citiesVisted[i] = 0;
			ants[ant].pathTaken[i] = -1;
		}
		if (cityTo == numberOfCities)
		{
			cityTo = 0;
		}
		ants[ant].currentCity = cityTo++;
		ants[ant].pathNumhber = 1;
		ants[ant].pathTaken[0] = ants[ant].currentCity;
		ants[ant].citiesVisted[ants[ant].currentCity] = 1;
	}

	tourAverage = tourTotal / (numberOfAnts*0.8);
	cout << "tour average " << tourAverage << endl;
	

}


double antInformation(int cityFrom, int cityTo)
{
	double info = ((pow(pheromones[cityFrom][cityTo], ALPHA)* pow((1.0 / distance[cityFrom][cityTo]), BETA)));
	//cout << info << endl;
	return info;

}


int selectNextCity(int ant)
{
	int citiesChecked = 0;
	int cityFrom, cityTo;
	double equationDenominator = 0.0;
	cityFrom = ants[ant].currentCity;
	for (cityTo = 0; cityTo < numberOfCities; cityTo++)
	{
		//cout << endl;
		//cout << "city from" << cityFrom << endl;
		//cout << "city to" << cityTo << endl;
		if (ants[ant].citiesVisted[cityTo] == 0)
		{
			equationDenominator += antInformation(cityFrom, cityTo);
			//cout << antInformation(cityFrom, cityTo) << endl;
			//cout << equationDenominator << endl;
		}
	}
	assert(equationDenominator != 0.0);
	do {
		double information;
		cityTo++;
		if (cityTo >= numberOfCities)
		{
			cityTo = 0;
		}
		if (ants[ant].citiesVisted[cityTo] == 0)
		{
			information = antInformation(cityFrom, cityTo) / equationDenominator;
			double x = ((double)rand() / RAND_MAX);
			if (x < information)
			{
				break;
			}
		}
		/*citiesChecked++;
		//cout << citiesChecked << endl;
		if (citiesChecked == numberOfCities)
		{
			break;
		}*/
	} while (1);
	return cityTo;
}


int antColonySimulation()
{
	int moving = 0;

	for (int i = 0; i < numberOfCities; i++)
	{
		if (ants[i].pathNumhber < numberOfCities)
		{
			//cout << ants[i].currentCity << ",";
			ants[i].nextCity = selectNextCity(i);
			//cout << ants[i].nextCity << ",";
			ants[i].citiesVisted[ants[i].nextCity] = 1;
			ants[i].pathTaken[ants[i].pathNumhber++] = ants[i].nextCity;
			//cout << ants[i].currentCity << ",";
			//cout << ants[i].nextCity << ",";
			ants[i].pathDistance += distance[ants[i].currentCity][ants[i].nextCity];

			if (ants[i].pathNumhber == numberOfCities)
			{
				ants[i].pathDistance += distance[ants[i].pathTaken[numberOfCities - 1]][ants[i].pathTaken[0]];
			}
			ants[i].currentCity = ants[i].nextCity;
			moving++;
			//cout << moving << endl;
		}

	}
	//cout << moving << endl;
	return moving;
}


void updateTrails()
{
	int cityFrom, cityTo, ant;

	for (cityFrom = 0; cityFrom < numberOfCities; cityFrom++)
	{
		for (cityTo = 0; cityTo < numberOfCities; cityTo++)
		{
			if (cityFrom != cityTo)
			{
				pheromones[cityFrom][cityTo] *= (1.0 - evaporationRateFactor);
				if (pheromones[cityFrom][cityTo] < 0.0)
				{
					pheromones[cityFrom][cityTo] = initial_Pheromones;
				}
			}
		}
	}



	//assign new pheromones
	for (ant = 0; ant < numberOfAnts; ant++)
	{
		for (int i = 0; i < numberOfCities; i++)
		{
			if (i < numberOfCities - 1)
			{
				cityFrom = ants[ant].pathTaken[i];
				cityTo = ants[ant].pathTaken[i + 1];
			}
			else
			{
				cityFrom = ants[ant].pathTaken[i];
				cityTo = ants[ant].pathTaken[0];
			}

			pheromones[cityFrom][cityTo] += (QVal / ants[ant].pathDistance);
			pheromones[cityTo][cityFrom] = pheromones[cityFrom][cityTo];
		}
	}

	for (cityFrom = 0; cityFrom < numberOfCities; cityFrom++)
	{
		for (cityTo = 0; cityTo < numberOfCities; cityTo++)
		{
			pheromones[cityFrom][cityTo] *= evaporationRateFactor;
		}
	}
}


void antColonyAlgorithm()
{
	cout << endl;
	cout << endl;
	cout << "ANT COLONY ALGORITHM" << endl;
	cout << endl;


	AntsType bestAnt;
	int cumulativePathTotal = 0;
	while (cumulativePathTotal++ < maxCombinations)
	{
		if (antColonySimulation() == 0)
		{
			updateTrails();

			if (cumulativePathTotal != maxCombinations)
			{
				restartAllAnts();
				//cout << count << endl;
				//count++;
			}
			cout << "current number of combinations tried " << cumulativePathTotal << endl;
			cout << "Current best is" << "(" << best << ")" << endl;
		}
	}
	bestAnt = ants[bestPathNumber];
	for (int i = 0; i < numberOfCities; i++)
	{
		cout << bestAnt.pathTaken[i] << ",";
	}
	cout << endl;
	//cout << "best path was = " << bestPathNumber << endl;
}


void geneticAlgorithm()
{
	cout << endl;
	cout << endl;
	cout << "Genetic ALGORITHM" << endl;
	cout << endl;
	genes.clear();
	//int count = 0;
	generateGenes();
	for (int i = 0; i < numberOfGenes; i++)
	{
		genes.push_back(startingGenes[i]);
	}



	calculateRouteDistance();
	for (int generation = 0; generation < numberOfGenerations-1; generation++)
	{
		selectGenes(genes, fitness);
		calculateRouteDistance();
		//count += 1;
		if (bestDistance < record)
		{
			record = bestDistance;
			//count = 0;
			cout << "new record" << endl;
		}

	}
}



void hybridAlgorithm()
{
	cout << endl;
	cout << endl;
	cout << "HYBRID ALGORITHM" << endl;
	cout << endl;
	genes.clear();
	best = (double)maxRouteLength;
	bestDistance = (double)maxRouteLength;
	record = (double)maxRouteLength;
	bestPathNumber = NULL;
	
	Hybrid = true;
	//int count = 0;
	if (antColonySimulation() == 0)
	{
		updateTrails();
		restartAllAnts();
		//cout << count << endl;
		//count++;
	}
	cout << "Current best is" << "(" << best << ")" << endl;
	/*for (int i = 0; i < numberOfGenes; i++)
	{
		for (int j = 0; j < numberOfCities; j++)
		{
			cout << genes[i].gene[j] << ", ";
		}
		cout << endl;
	}*/
	calculateRouteDistance();
	for (int generation = 0; generation < numberOfGenerations-1; generation++)
	{
		selectGenes(genes, fitness);
		calculateRouteDistance();
		//count += 1;
		if (bestDistance < record)
		{
			record = bestDistance;
			//count = 0;
			cout << "new record" << endl;
		}

	}
}



int main()
{
	//int count = 0;

	srand(time(NULL));

	initiate();

/*
	auto startGenerateGene = std::chrono::high_resolution_clock::now();
	
	auto endGenerateGene = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeTakenGenerate = endGenerateGene - startGenerateGene;
	cout << endl;
	cout << "time taken for Generate Genes: " << timeTakenGenerate.count() << endl;
	cout << endl;
	*/

	cout << "total number of combinations that will be tried = " << maxCombinations << endl;;
	double totalPossibleCombinations = 1;

	for (int i = 1; i < numberOfCities; i++)
	{
		totalPossibleCombinations = totalPossibleCombinations * i;
	}
	totalPossibleCombinations = totalPossibleCombinations / 2;
	cout << "total possible orders = " << totalPossibleCombinations << endl;






	auto startAnt = std::chrono::high_resolution_clock::now();
	antColonyAlgorithm();
	auto finishAnt = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeTakenAnt = finishAnt - startAnt;
	cout << endl;
	cout << "time taken for Ant: " << timeTakenAnt.count() << endl;
	cout << endl;





	auto startGenetic = std::chrono::high_resolution_clock::now();
	geneticAlgorithm();
	auto finishGenetic = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeTakenGenetic = finishGenetic - startGenetic;
	cout << endl;
	cout << "time taken for Genetic Algorithm: " << timeTakenGenetic.count() << endl;
	cout << endl;






	auto startHybrid = std::chrono::high_resolution_clock::now();
	hybridAlgorithm();
	auto finishHybrid = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeTakenHybrid = finishHybrid - startHybrid;
	cout << endl;
	cout << "time taken for Hybrid: " << timeTakenHybrid.count() << endl;
	cout << endl;
	
	return 0;
}

