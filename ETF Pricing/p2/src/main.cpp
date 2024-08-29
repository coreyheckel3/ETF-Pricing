/*
Problem 1: (50 points)

    In the data folder, there is a data.csv containing 6 month daily pricing data for 24 ETF stocks, the pricing data follows
    the chronological order (from earliest to the latest).

    1. (10 points) read in the ETF pricing information from the csv file and store the data in a two-dimensional vector.
    2. (10 points) Calculate the daily returns for each ETF on each date.
    2. (10 points) Calculate the mean daily returns for each ETF.
    3. (10 points) Calculate the daily volatility for each ETF.
    2. (10 points) calculate the covariance matrix for those stock returns and print out the covariance matrix.

    Hint 0: I provide some function prototypes for you, feel free to use it and feel free to create your own functions.

    Hint 1: If the shape of stock pricing matrix is with the shape (m,n), then the stock return matrix's shape is (m-1,n).

    Hint 2: For covariance matrix M: M(i,j) = Cov(stockReturn_i, stockReturn_j). For Cov(x,y) = E(xy) - E(x)E(y).

    Hint 3: You may need to transpose the return matrix or just define the return matrix shape to be (n, m-1) when you are
        calculating the returns, so that every row of the matrix contains all the daily returns for a single stock. This will
        make it easy when you are calculating the covariance matrix.

    Hint 4: You do not need to use any of the previous hints to successfully complete this task.
*/



#include "std_lib_facilities.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

vector<vector<double>> readData(string fileName) {
    ifstream file(fileName);
    if (!file.is_open()) {
        cout << "Failed to open file: " << fileName << endl;
        exit(1);
    }

    vector<vector<double>> data;
    string line;
    while (getline(file, line)) {
        vector<double> row;
        stringstream ss(line);
        string cell;

        while (getline(ss, cell, ',')) {
            row.push_back(stod(cell));
        }

        data.push_back(row);
    }

    file.close();
    return data;
}
// Calculate the daily returns of all ETFs
vector<vector<double>> stockReturn(vector<vector<double>>& stockPrice){
    vector<vector<double>> sReturn(stockPrice[0].size(), vector<double>(stockPrice.size()-1));

    for (int i=0; i<stockPrice.size()-1; i++){
        for (int j=0; j<stockPrice[0].size(); j++){
            sReturn[j][i] = (stockPrice[i+1][j] - stockPrice[i][j]) / stockPrice[i][j];
        }
    }

    return sReturn;
}

// Calculate the mean of all returns
vector<double> meanReturn(vector<vector<double>>& stockReturn){
    vector<double> mean(stockReturn.size());

    for (int i=0; i<stockReturn.size(); i++){
        double sum = 0;
        for (int j=0; j<stockReturn[0].size(); j++){
            sum += stockReturn[i][j];
        }
        mean[i] = sum / stockReturn[0].size();
    }

    return mean;
}

// Calculate the volatility of all ETFs
vector<double> dailyVolatility(vector<vector<double>>& stockReturn){
    vector<double> volatility(stockReturn.size());

    for (int i=0; i<stockReturn.size(); i++){
        double sum = 0;
        double mean = 0;
        for (int j=0; j<stockReturn[0].size(); j++){
            mean += stockReturn[i][j];
        }
        mean /= stockReturn[0].size();
        for (int j=0; j<stockReturn[0].size(); j++){
            sum += pow(stockReturn[i][j]-mean, 2);
        }
        volatility[i] = sqrt(sum / (stockReturn[0].size()-1));
    }

    return volatility;
}

// Calculate covariance between two stocks
double Cov(vector<double>& stockReturn1, vector<double>& stockReturn2){
    double mean1 = 0, mean2 = 0;
    for (int i=0; i<stockReturn1.size(); i++){
        mean1 += stockReturn1[i];
        mean2 += stockReturn2[i];
    }
    mean1 /= stockReturn1.size();
    mean2 /= stockReturn1.size();

    double cov = 0;
    for (int i=0; i<stockReturn1.size(); i++){
        cov += (stockReturn1[i]-mean1) * (stockReturn2[i]-mean2);
    }

    cov /= stockReturn1.size()-1;

    return cov;
}

// Calculate covariance matrix of the stocks
vector<vector<double>> covarianceMatrix(vector<vector<double>>& stockReturn){
    vector<vector<double>> cMatrix(stockReturn.size(), vector<double>(stockReturn.size()));

    for (int i=0; i<stockReturn.size(); i++){
        for (int j=0; j<stockReturn.size(); j++){
            cMatrix[i][j] = Cov(stockReturn[i], stockReturn[j]);
        }
    }

    return cMatrix;
}

// helper functions use to print out vectors and matrices
void printVector(vector<double>& vect){
    for (double v:vect){
        cout << v << " ";
    }
    cout << endl;
}

void printMatrix(vector<vector<double>>& vect) {
    for (vector<double> v: vect) {
        for (double vv: v) {
            cout << vv << " ";
        }
        cout << endl;
    }
}
    int main() {
        ifstream inFile;
        string fileName = "C:/Users/wheel/Downloads/midtermExam/midtermExam/p2/data/data.csv";
        vector<vector<double>> data = readData(fileName);


        for (const auto& row : data) {
            for (double val : row) {
                cout << val << ",";
            }
            cout << endl;
        }
        const int N = 24;

        // Read in tickers
        string tickers[N];
        string line;
        getline(inFile, line); // read the first line which contains the tickers
        istringstream ss(line);
        string token;
        int i = 0;
        while (getline(ss, token, ',')) {
            tickers[i] = token;
            i++;
        }

// Read in stock prices
        vector<vector<double>> stockPrice;
        while (getline(inFile, line)) {
            vector<double> prices;
            istringstream ss(line);
            while (getline(ss, token, ',')) {
                prices.push_back(stod(token));
            }
            stockPrice.push_back(prices);
        }

// calculate the stock returns
        vector<vector<double>> stockReturn(N, vector<double>(stockPrice.size() - 1));
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < stockPrice.size() - 1; j++) {
                double r = (stockPrice[j + 1][i] - stockPrice[j][i]) / stockPrice[j][i];
                stockReturn[i][j] = r;
            }
        }

// calculate and print out the mean returns
        vector<double> meanReturns;
        for (int i = 0; i < N; i++) {
            double mean = accumulate(stockReturn[i].begin(), stockReturn[i].end(), 0.0) / stockReturn[i].size();
            meanReturns.push_back(mean);
        }
        cout << "Mean Returns: " << endl;
        printVector(meanReturns);

// calculate and print out the stocks' volatility
        vector<double> volatility;
        for (int i = 0; i < N; i++) {
            double sq_sum = inner_product(stockReturn[i].begin(), stockReturn[i].end(), stockReturn[i].begin(), 0.0);
            double var = sq_sum / stockReturn[i].size() - meanReturns[i] * meanReturns[i];
            volatility.push_back(sqrt(var));
        }
        cout << "Volatility: " << endl;
        printVector(volatility);

// calculate and print out the covariance matrix
        vector<vector<double>> covariance = covarianceMatrix(stockReturn);
        cout << "Covariance Matrix: " << endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                cout << covariance[i][j] << " ";
            }
            cout << endl;
        }

        inFile.close();
        return 0;
    }