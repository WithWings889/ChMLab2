#include <iostream>
#include <vector>
#include <string>

using namespace std;

void showProblem(vector<vector<double>> a, vector<double> y, int n)
{

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << a[i][j] << "*x" << j;
			if (j < n - 1)
				cout << " + ";
		}
		cout << " = " << y[i] << endl;
	}
}
vector<vector<double>> matrixMultip(vector<vector<double>> &a, vector<vector<double>> &b )
{
	vector<vector<double>> c(a.size(), vector<double>(b[0].size()));
	for (int i = 0; i < a.size(); i++)
	{
		//c[i] = vector<double>(a[0].size());
		for (int j = 0; j < b[0].size(); j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < a[0].size(); k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
	return c;
}

vector<double> matrixMultip(vector<vector<double>>& a, vector<double>& b)
{
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); i++)
	{
		c[i] = 0;
		//c[i] = vector<double>(a[0].size());
		for (int j = 0; j < b.size(); j++)
		{
			c[i] += a[i][j] * b[j];
		}
	}
	return c;
}
void inicialization(vector<vector<double>> &a, vector<double> &y, int n)
{

	vector<double> tempVector;
	for (int i = 0; i < n; i++)
	{
		double temp;
		for (int j = 0; j < n; j++)
		{
			cout << "a[" << i << "][" << j << "]= ";
			cin >> temp;
			tempVector.push_back(temp);
		}
		a.push_back(tempVector);
		tempVector.clear();
	}
	for (int i = 0; i < n; i++)
	{
		double temp;
		cout << "y[" << i << "]= ";
		cin >> temp;
		y.push_back(temp);
	}
}

vector<double> Gauss(vector<vector<double>> a, vector<double> y, int n, double &detA)
{
	double max;
	vector<double> x(n);
	int k, index;
	const double eps = 0.00001;
	k = 0;
	while (k < n)
	{
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		try
		{
			if (max < eps)
			{
				string error = "Рішення неможливо отримати через нульовий стовпчик " + to_string(index) + " матриці A\n";
				throw error;
			}
		}
		catch (string error)
		{
			cout << error;
		}
		detA *= max;
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;

		vector<vector<double>> m(n, vector<double>(n));
		for (int i = 0; i < n; ++i)
		{
			m[i][i] = 1;
		}

		for (int i = k; i < n; i++)
		{
			if (abs(a[i][k]) < eps)
				continue;

			if (i == k)
				m[i][k] = 1 / a[k][k];
			else
				m[i][k] = (-1) * a[i][k] / a[k][k];
		}

			a = matrixMultip(m, a);
			y = matrixMultip(m, y);
		
		k++;
	}
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}
void transponation(vector<vector<double>>& a, int n)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			double temp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = temp;
		}
	}
}
double norm(vector<vector<double>> a, int n)
{
	double max = 0, sum = 0;
	for (int i = 0; i < n; ++i)
		max += abs(a[0][i]);
	for (int j = 0; j < n; ++j)
	{
		for (int i = 1; i < n; ++i)
		{
			sum += abs(a[i][j]);
		}
		if (max < sum)
		{
			max = sum;
		}
		sum = 0;
	}

	return max;
}

double conditionNumber(vector<vector<double>> a, int n)
{
	//Обернена
	double detA = 1;
	vector<vector<double>> b;
	vector<double> e(n);
	for (int i = 0; i < n; ++i)
	{
		e[i] = 1;
		b.push_back(Gauss(a, e, n, detA));
		e[i] = 0;
	}
	transponation(b, n);

	return norm(a, n) * norm(b, n);
}

double findQ(vector<vector<double>> a, int n)
{
	double max = 0, sum = 0;
	for (int i = 1; i < n; ++i)
	{
		max += abs(a[0][i]);
	}
	max /= abs(a[0][0]);
	for (int j = 1; j < n; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			if (i == j) continue;
			sum += abs(a[j][i]);
		}
		sum /= abs(a[j][j]);
		if (max < sum)
		{
			max = sum;
		}
		sum = 0;
	}
	return max;
}
double normaVector(vector<double> x1, vector<double> x2)
{
	if (x1.size() != x2.size())
		return -1;
	double max = abs(x1[0] - x2[0]);
	for (int i = 1; i < x1.size(); ++i)
	{
		if (max < abs(x1[i] - x2[i]))
			max = abs(x1[i] - x2[i]);
	}
	return max;
}

int apriorYakobi(vector<vector<double>> a, vector<double> y, int n, vector<double> x0, double e)
{
	double sum = 0, q = findQ(a, n);

	//find x1
	vector<double> x1;
	for (int k = 0; k < n; ++k)
	{
		for (int i = 0; i < n; ++i)
		{
			if (i == k) continue;
			sum -= a[i][k] * x0[i];
		}
		sum = (sum + y[k]) / a[k][k];
		x1.push_back(sum);
	}

	double ans = log(normaVector(x0, x1) / ((1 - q) * e)) / log(1 / q);
	return trunc(ans) + 1;
}



vector<double> Yakobi(vector<vector<double>> a, vector<double> y, double e, int n, int& apostarior)
{
	vector<double> x(n);
	vector<double> xNew;
	double sum = 0;
	while(true)
	{
		for (int k = 0; k < n; ++k)
		{
			for (int i = 0; i < n; ++i)
			{
				if (i == k) continue;
				sum -= a[i][k] * x[i];
			}
			sum = (sum + y[k]) / a[k][k];
			xNew.push_back(sum);
			sum = 0;
		}
		++apostarior;
		if (normaVector(x, xNew) < e)
			break;
		x = xNew;
		xNew.clear();
	} 
	return x;
}

	void showAnsGauss(vector<vector<double>> a, vector<double> y, int n)
	{
		double detA = 1;
		vector<double> x;
		x = Gauss(a, y, n, detA);
		cout << "-----------------Gauss--------------------\n";
		for (int i = 0; i < n; i++)
			cout << "x[" << i << "]=" << x[i] << endl;
		cout << "cond(A) is " << conditionNumber(a, n) << endl;
		cout << "det A is " << detA << endl;
		cout << "------------------------------------------\n\n";
	}
	bool YakobiIsAvailable(vector<vector<double>> a)
	{
		double sum = 0;
		for (int i = 0; i < a.size(); ++i)
		{
			for (int j = 0; j < a.size(); ++j)
			{
				if (j == i)
					continue;
				else
				{
					sum += a[j][i];
				}
			}
			if (a[i][i] < sum)
				return false;
			sum = 0;
		}
		return true;
	}
	void showAnsYakobi(vector<vector<double>> a, vector<double> y, int n, double e)
	{
		vector<double> x;
		vector<double> x0(n);
		int apostarior = 0;
		if (YakobiIsAvailable(a))
		{
			x = Yakobi(a, y, e, n, apostarior);
			cout << "-----------------Yakobi--------------------\n";
			for (int i = 0; i < n; i++)
				cout << "x[" << i << "]=" << x[i] << endl;
			cout << "Apostarior is " << apostarior << endl;;
			cout << "Aprior is " << apriorYakobi(a, y, n, x0, e) << endl;
			cout << "q = " << findQ(a, n) << endl;
			cout << "--------------------------------------------\n";
		}
		else
			cout << "Yakobi method is not available" << endl;
		cin.get(); cin.get();
	}

	int main()
	{
		int n;
		double e = 1e-2;
		system("chcp 1251");
		system("cls");
		cout << "Ведіть кількість рівнянь: ";
		cin >> n;
		vector<vector<double>> a;
		vector<double>y;
		inicialization(a, y, n);

		showProblem(a, y, n);
		
		showAnsGauss(a, y, n);
		showAnsYakobi(a, y, n, e);


		return 0;
	}