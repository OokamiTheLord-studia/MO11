#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

using namespace std;

//Definicja struktury node
struct node {
	double t;
	double x;
	double value;
};


///Zadanie 6
class DU
{

public:
	vector<double> d;
	vector<double> u;

	const int siz;
	DU(int size) :siz(size)
	{
		d.resize(siz);
		u.resize(siz - 1);
	}

	double at(unsigned int i, unsigned int j)
	{

		if (i < siz && j < siz)
		{
			if (i == j) return d[j];
			if (j == i + 1) return u[j];
			return 0;
		}
		else
		{
			throw invalid_argument("numery wierszy i kolumn musza byc z przedzialu <0,rozmiar)");
		}
	}

	/*void zapisz(unsigned int i, unsigned int j, double w)
	{
		if (i < siz && j < siz)
		{

			if (i == j) { d[j] = w; }
			else if (j == i + 1) { u[j] = w; }
			else throw invalid_argument("To nie jest lokalizacja mozliwa do zapisu");

		}
		else
		{
			throw invalid_argument("numery wierszy i kolumn musza byc z przedzialu <0,rozmiar)");
		}
	}*/

	friend class Tridiagonal;
	friend void solver(Tridiagonal& tri, DU& du, vector<double>& b, vector<double>& x);
};


class Tridiagonal
{
	vector<double> l;
	vector<double> d;
	vector<double> u;


public:
	const int siz;
	Tridiagonal(vector<double> low, vector<double> diag, vector<double> upp) : l(low), d(diag), u(upp), siz(diag.size())
	{
		if ((low.size() + 1 != diag.size()) || (upp.size() + 1 != diag.size())) throw invalid_argument("zla dlugosc macierzy low badz upp");
	}

	double at(unsigned int i, unsigned int j)
	{

		if (i < siz && j < siz)
		{
			if (i == j) return d[j];
			if (j == i + 1) return u[j];
			if (j == i - 1) return l[j];
			return 0;
		}
		else
		{
			throw invalid_argument("numery wierszy i kolumn musza byc z przedzialu <0,rozmiar)");
		}
	}

	void ThomasMatrix(DU& du)
	{
		if (du.siz != siz) throw invalid_argument("macierz LU musi miec ten sam rozmiar");
		du.u = u;
		du.d[0] = d[0];
		for (int i = 1; i < siz; i++)
		{
			du.d[i] = d[i] - ((l[i - 1] * u[i - 1]) / du.d[i - 1]);
		}

	}

	friend void solver(Tridiagonal& tri, DU& du, vector<double>& b, vector<double>& x);
	friend void node_solver(Tridiagonal& tri, DU& du, vector<double>& b, vector<node>* x);
};

void solver(Tridiagonal& tri, DU& du, vector<double>& b, vector<double>& x)
{
	const int locsize = tri.siz;

	if ((locsize != du.siz) || (locsize != b.size())) throw invalid_argument("wektory i macierze musza miec ten sam rozmiar");
	x.resize(locsize);

	vector<double> r;
	r.resize(locsize);

	r[0] = b[0];
	for (int i = 1; i < locsize; i++)
	{
		r[i] = b[i] - ((tri.l[i - 1] * r[i - 1]) / du.d[i - 1]);
	}

	x[locsize - 1] = r[locsize - 1] / tri.d[locsize - 1];
	for (int i = locsize - 2; i >= 0; i--)
	{
		x[i] = (r[i] - (tri.u[i] * x[i + 1])) / tri.d[i];
	}

}


void node_solver(Tridiagonal& tri, DU& du, vector<double>& b, vector<node>* x)
{
	const int locsize = tri.siz;

	if ((locsize != du.siz) || (locsize != b.size())) throw invalid_argument("wektory i macierze musza miec ten sam rozmiar");
	//x.resize(locsize);

	vector<double> r;
	r.resize(locsize);

	r[0] = b[0];
	for (int i = 1; i < locsize; i++)
	{
		r[i] = b[i] - ((tri.l[i - 1] * r[i - 1]) / du.d[i - 1]);
	}

	//x[locsize - 1] = r[locsize - 1] / tri.d[locsize - 1];
	x->at(locsize -1).value = r[locsize - 1] / tri.d[locsize - 1];
	for (int i = locsize - 2; i >= 0; i--)
	{
		x->at(i).value = (r[i] - (tri.u[i] * x->at(i + 1).value)) / tri.d[i];
	}

}


///Zadanie 7
class Wector
{
private:
	std::vector<double> body;

public:
	unsigned int siz;

	Wector(int rozm) : siz(rozm)
	{
		body.resize(rozm);
	}

	Wector(initializer_list<double> ilist, unsigned int size) : siz(size)
	{
		if (ilist.size() != size) throw invalid_argument("rozmiar nie zgadza sie z iloscia elementow!");
		body = ilist;
	}


	void Add(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby dodawac wektory musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < result.siz; i++)
			result.zapisz(i) = V.at(i) + this->at(i);
	}

	void Sub(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby odejmowac wektory musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < result.siz; i++)
			result.zapisz(i) = this->at(i) - V.at(i);
	}

	void SubAbs(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby odejmowac wektory musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < result.siz; i++)
			result.zapisz(i) = abs(this->at(i) - V.at(i));
	}

	void Abs()
	{
		for (unsigned int i = 0; i < body.size(); i++)
			body[i] = abs(body[i]);
	}

	void Abs(Wector& V)
	{
		if (siz != V.siz) throw invalid_argument("Wektory musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < body.size(); i++)
			V.body[i] = abs(body[i]);
	}

	const double Max()
	{
		return *max_element(body.begin(), body.end());
	}

	const double& at(unsigned int i)
	{
		return body[i];
	}

	double& zapisz(unsigned int i)
	{
		return body[i];
	}
};

//Klasa Matrix jest tylko interfejsem
class FullMatrix;

class Matrix
{

protected:
	const double zero = 0;

public:
	const unsigned int siz;

	Matrix(int siz) :siz(siz) {}



	virtual const double& at(unsigned int i, unsigned int j) = 0;

	virtual void Multiply(Wector& V, Wector& result) = 0;
	virtual void Multiply(Matrix& M, FullMatrix& result) = 0;

	virtual void Add(Matrix& M, FullMatrix& result) = 0;

};



class FullMatrix :public Matrix
{
private:
	vector<double> body;

public:
	FullMatrix(int rozm) : Matrix(rozm)
	{
		body.resize(rozm * rozm);
	}

	FullMatrix(initializer_list<double> ilist, unsigned int size) : Matrix(size)
	{
		if ((int)sqrt(ilist.size()) != size) throw invalid_argument("rozmiar nie zgadza sie z iloscia elementow! Macierze moga byc tylko kwadratowe");
		body = ilist;
	}

	void Multiply(Matrix& M, FullMatrix& result)
	{
		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierze musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
			{
				double temp = 0;
				for (unsigned int z = 0; z < result.siz; z++)
					temp += this->at(j, z) * M.at(z, i);
				result.zapisz(j, i) = temp;
			}


	}

	void Multiply(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierz i wektory musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < V.siz; i++)
		{
			double temp = 0;
			for (unsigned int z = 0; z < V.siz; z++)
				temp += this->at(i, z) * V.at(z);
			result.zapisz(i) = temp;
		}


	}

	void Add(Matrix& M, FullMatrix& result)
	{

		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby dodawac macierze musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
				result.zapisz(j, i) = this->at(j, i) + M.at(j, i);
	}

	void Negate()
	{
		for (unsigned int i = 0; i < body.size(); i++)
			body[i] = -body[i];
	}

	void Negate(FullMatrix& result)
	{
		if ((this->siz != result.siz) || (this->siz != result.siz)) throw invalid_argument("Macierze musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < body.size(); i++)
			result.body[i] = -body[i];
	}


	const double& at(unsigned int i, unsigned int j)
	{
		return body[i * siz + j];
	}

	double& zapisz(unsigned int i, unsigned int j)
	{
		return body[i * siz + j];
	}


};

class TUpper :public Matrix
{
private:
	std::vector<double> body;

public:
	TUpper(int rozm) : Matrix(rozm)
	{
		int temp = (rozm - 1) * rozm / 2;
		body.reserve(temp);
		body.resize(temp);
	}



	TUpper(FullMatrix& m) : Matrix(m.siz)
	{
		int temp = (siz - 1) * siz / 2;
		body.reserve(temp);
		body.resize(temp);

		unsigned int ind = 0;
		for (unsigned int i = 1; i < siz; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				body[ind] = m.at(j, i);
				ind++;
			}
		}
	}

	void Multiply(Matrix& M, FullMatrix& result)
	{
		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierze musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = result.siz - 1; j < UINT32_MAX; j--)
			{
				double temp = 0;
				for (unsigned int z = result.siz - 1; z > j; z--)
					temp += this->at(j, z) * M.at(z, i);
				result.zapisz(j, i) = temp;
			}


	}

	void Multiply(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierz i wektory musza byc tych samych rozmiarow");

		for (unsigned int i = V.siz - 1; i < UINT32_MAX; i--)
		{
			double temp = 0;
			for (unsigned int z = V.siz - 1; z > i; z--)
				temp += this->at(i, z) * V.at(z);
			result.zapisz(i) = temp;
		}
	}

	void Add(Matrix& M, FullMatrix& result)
	{
		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby dodawac macierze musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
				result.zapisz(j, i) = this->at(j, i) + M.at(j, i);
	}

	const double& at(unsigned int i, unsigned int j)
	{
		if (i >= j)
		{
			return zero;
		}
		else
		{
			return body[j * (j - 1) / 2 + i];
		}
	}

	double& zapisz(int i, int j)
	{
		if (i >= j)
		{
			throw std::invalid_argument("Nie mozna pisac do elementow macierzy trojkatnej gornej powyzej przekatnej!!!");
		}
		else
		{
			return (body[j * (j - 1) / 2 + i]);
		}
	}
};


class Diagonal :public Matrix
{
private:
	vector<double> body;

public:
	Diagonal(unsigned int rozmiar) : Matrix(rozmiar)
	{
		body.reserve(siz);
		body.resize(siz);
	}

	Diagonal(FullMatrix& m) :Matrix(m.siz)
	{
		body.reserve(siz);
		body.resize(siz);

		for (unsigned int i = 0; i < siz; i++)
		{
			body[i] = m.at(i, i);
		}
	}

	void Multiply(Matrix& M, FullMatrix& result)
	{

		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierze musza byc tych samych rozmiarow");


		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
			{
				result.zapisz(j, i) = body[j] * M.at(j, i);
			}


	}

	void Multiply(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierz i wektory musza byc tych samych rozmiarow");

		for (unsigned int i = 0; i < V.siz; i++)
			result.zapisz(i) = body[i] * V.at(i);
	}

	void Add(Matrix& M, FullMatrix& result)
	{
		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby dodawac macierze musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
				result.zapisz(j, i) = this->at(j, i) + M.at(j, i);
	}

	void Invert()
	{
		for (unsigned int i = 0; i < body.size(); i++)
			body[i] = 1 / body[i];
	}

	void Invert(Diagonal& D)
	{
		if (this->siz != D.siz) throw invalid_argument("Macierze musza byc rownych rozmiarow");
		for (unsigned int i = 0; i < body.size(); i++)
			D.body[i] = 1 / body[i];
	}

	const double& at(unsigned int i, unsigned int j)
	{
		if (i == j)
		{
			return body[i];
		}
		else
		{
			return zero;
		}
	}

	double& zapisz(int i, int j)
	{
		if (i == j)
		{
			return (body[i]);
		}
		else
		{
			throw std::invalid_argument("Nie mozna pisac do elementow macierzy trojkatnej gornej powyzej przekatnej!!!");

		}
	}

};


class TLower :public Matrix
{
private:
	std::vector<double> body;

public:
	TLower(int rozm) : Matrix(rozm)
	{
		int temp = (rozm - 1) * rozm / 2;
		body.reserve(temp);
		body.resize(temp);
	}



	TLower(FullMatrix& m) : Matrix(m.siz)
	{
		int temp = (siz - 1) * siz / 2;
		body.reserve(temp);
		body.resize(temp);

		unsigned int ind = 0;
		for (unsigned int j = 1; j < siz; j++)
		{
			for (unsigned int i = 0; i < j; i++)
			{
				body[ind] = m.at(j, i);
				ind++;
			}
		}
	}


	void Multiply(Matrix& M, FullMatrix& result)
	{
		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierze musza byc tych samych rozmiarow");



		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
			{
				double temp = 0;
				for (unsigned int z = 0; z < j; z++)
					temp += this->at(j, z) * M.at(z, i);
				result.zapisz(j, i) = temp;
			}


	}

	void Multiply(Wector& V, Wector& result)
	{
		if ((this->siz != V.siz) || (this->siz != result.siz)) throw invalid_argument("Aby mnozyc macierz i wektory musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < V.siz; i++)
		{
			double temp = 0;
			for (unsigned int z = 0; z < i; z++)
				temp += this->at(i, z) * V.at(z);
			result.zapisz(i) = temp;
		}
	}

	void Add(Matrix& M, FullMatrix& result)
	{
		if ((this->siz != M.siz) || (this->siz != result.siz)) throw invalid_argument("Aby dodawac macierze musza byc tych samych rozmiarow");
		for (unsigned int i = 0; i < result.siz; i++)
			for (unsigned int j = 0; j < result.siz; j++)
				result.zapisz(j, i) = this->at(j, i) + M.at(j, i);
	}



	const double& at(unsigned int i, unsigned int j)
	{
		if (i <= j)
		{
			return zero;
		}
		else
		{
			return body.at(i * (i - 1) / 2 + j);
		}
	}

	double& zapisz(int i, int j)
	{
		if (i < j)
		{
			throw std::invalid_argument("Nie mozna pisac do elementow macierzy trojkatnej gornej powyzej przekatnej!!!");
		}
		else
		{
			return (body[i * (i - 1) / 2 + j]);
		}
	}
};


void Jakobi(FullMatrix& A, Wector& x0, Wector& b, double etol, double ftol, unsigned int iter, Wector& solved)
{
	Diagonal Dinv(A);
	Dinv.Invert();
	TLower L(A);
	TUpper U(A);
	Wector C(Dinv.siz);
	Dinv.Multiply(b, C);
	FullMatrix LU(A.siz);
	L.Add(U, LU);
	FullMatrix M(A.siz);
	Dinv.Multiply(LU, M);
	M.Negate();

	Wector x(x0.siz);
	Wector xx(x0.siz);

	M.Multiply(x0, xx);
	xx.Add(C, xx);

	Wector E(xx.siz);
	xx.SubAbs(x0, E);
	Wector F(solved.siz);
	A.Multiply(xx, F);
	F.SubAbs(b, F);
	unsigned int n = 0;

	printf_s("ITERATION %d\n", n);
	for (unsigned int i = 0; i < xx.siz; i++)
		printf("%9.6e ", xx.at(i));
	printf("\nERROR: %9.6e RESIDUUM %9.6e\n\n", E.Max(), F.Max());


	Wector* x1 = &x;
	Wector* x2 = &xx;

	while ((n < iter) && (E.Max() > etol) && (F.Max() > ftol))
	{
		swap(x1, x2);
		M.Multiply(*x1, *x2);
		(*x2).Add(C, *x2);

		(*x2).SubAbs(*x1, E);
		A.Multiply(*x2, F);
		F.SubAbs(b, F);
		n++;
		printf_s("ITERATION %d\n", n);
		for (unsigned int i = 0; i < xx.siz; i++)
			printf("%9.6e ", xx.at(i));
		printf("\nERROR: %9.6e RESIDUUM %9.6e\n\n", E.Max(), F.Max());
	}

	solved = (*x2);

}

void Gauss(FullMatrix& A, Wector& x0, Wector& b, double etol, double ftol, unsigned int iter, Wector& solved)
{
	TUpper U(A);
	FullMatrix LU(A.siz);
	TLower L(A);
	Diagonal Dinv(A);
	Dinv.Invert();
	L.Add(U, LU);

	Wector x(x0.siz);
	x = x0;
	unsigned int n = -1;
	Wector x1(x0.siz);
	x1 = x0;
	Wector E(x0.siz);
	Wector F(x0.siz);
	do
	{
		for (unsigned int i = 0; i < x.siz; i++)
		{
			double temp = 0;
			for (unsigned int j = 0; j < x.siz; j++)
				temp += LU.at(i, j) * x.at(j);
			x.zapisz(i) = Dinv.at(i, i) * (b.at(i) - temp);
		}
		x.SubAbs(x1, E);
		A.Multiply(x, F);
		F.SubAbs(b, F);
		n++;

		printf_s("ITERATION %d\n", n + 1);
		for (unsigned int i = 0; i < x.siz; i++)
			printf("%9.6e ", x.at(i));
		printf("\nERROR: %9.6e RESIDUUM %9.6e\n\n", E.Max(), F.Max());
		x1 = x;
	} while ((n < iter) && (E.Max() > etol) && (F.Max() > ftol));

	solved = x;
}

void SOR(FullMatrix& A, Wector& x0, Wector& b, double omega, double etol, double ftol, unsigned int iter, Wector& solved)
{
	TUpper U(A);
	FullMatrix LU(A.siz);
	TLower L(A);
	Diagonal Dinv(A);
	Dinv.Invert();
	L.Add(U, LU);
	double nomega = 1 - omega;

	Wector x(x0.siz);
	x = x0;
	unsigned int n = -1;
	Wector x1(x0.siz);
	x1 = x0;
	Wector E(x0.siz);
	Wector F(x0.siz);
	do
	{
		for (unsigned int i = 0; i < x.siz; i++)
		{
			double temp = 0;
			for (unsigned int j = 0; j < x.siz; j++)
				temp += LU.at(i, j) * x.at(j);
			x.zapisz(i) = Dinv.at(i, i) * omega * (b.at(i) - temp) + nomega * x.at(i);
		}
		x.SubAbs(x1, E);
		A.Multiply(x, F);
		F.SubAbs(b, F);
		n++;
		printf_s("ITERATION %d\n", n + 1);
		for (unsigned int i = 0; i < x.siz; i++)
			printf("%9.6e ", x.at(i));
		printf("\nERROR: %9.6e RESIDUUM %9.6e\n\n", E.Max(), F.Max());
		x1 = x;
	} while ((n < iter) && (E.Max() > etol) && (F.Max() > ftol));

	solved = x;
}

//Zadanie 11


//a musi być większe od 6
const int a = 7;
const int tmax = 1;
const int D = 1;
const double b = 0.1;

double rozw_analityczne(double x, double t)
{
	return 0.5 * exp(((D * t) / (b * b)) - (x / b)) * erfc((((2 * D * t) / (b)-x)) / (2 * sqrt(D * t)));
}

double warunek_poczatkowy(double x)
{
	return x < 0 ? 0 : exp(-x / b);
}

//struct node {
//	double t;
//	double x;
//	double value;
//};

class net
{
private:
	vector<vector<node>*> body;

public:
	net(int hcount, double h, int deltatcount, double deltat)
	{
		body.reserve(deltatcount);

		double current_t = 0.;

		for (int i = 0; i < deltatcount; i++)
		{
			vector<node>* temporary = new vector<node>;
			temporary->reserve(hcount + 2);

			double current_h = -a;

			node firstnode = { current_t, current_h - h, NAN };
			temporary->push_back(firstnode);

			for (int j = 0; j < hcount; j++)
			{

				node tempnode = { current_t, current_h, NAN };

				temporary->push_back(tempnode);

				current_h += h;
			}

			node lastnode = { current_t, current_h + h, NAN };
			temporary->push_back(lastnode);

			body.push_back(temporary);

			current_t += deltat;
		}
	}

	~net()
	{
		for (auto it = body.begin(); it != body.end(); it++)
		{
			delete* it;
		}
	}

	node* const getnode(int x, int y)
	{
		return &(body[x]->at(y));
	}

	size_t getxsize()
	{
		return body.size();
	}

	size_t getysize()
	{
		return body[0]->size();
	}

	friend net* crankanicolson_thomasa(double h, double deltat);
};



net* bezposrednia_eulera(double h, double deltat)
{
	const double lambda = D * (deltat / (h * h));
	assert(0.5 >= lambda);
	assert(0.43 > lambda);
	assert(0.37 < lambda);


	int hcount = (2. * a) / h;
	int deltatcount = tmax / deltat;

	//inicjalizacja
	net* localnet = new net(hcount, h, deltatcount, deltat);
	for (int i = 0; i < localnet->getysize(); i++)
	{
		auto tnode = localnet->getnode(0, i);
		tnode->value = warunek_poczatkowy(tnode->x);
	}

	//rozwiązanie
	for (int i = 1; i < localnet->getxsize(); i++)
	{
		for (int j = 1; j < localnet->getysize() - 1; j++)
		{
			auto node1 = localnet->getnode(i - 1, j - 1);
			auto node2 = localnet->getnode(i - 1, j);
			auto node3 = localnet->getnode(i - 1, j + 1);

			auto tnode = localnet->getnode(i, j);
			tnode->value = (lambda * node1->value) + ((1 - (2 * lambda)) * node2->value) + (lambda * node3->value);
		}

		//warunki brzegowe
		localnet->getnode(i, 0)->value = 0;
		localnet->getnode(i, localnet->getysize() - 1)->value = 0;
	}

	return localnet;

}

void dumpnet(net* paramnet, const char* filename)
{
	FILE* plik;
	fopen_s(&plik, filename, "w");

	fprintf_s(plik, ",");
	for (int i = 0; i < paramnet->getysize(); i++)
	{
		fprintf_s(plik, "%.16lf,", paramnet->getnode(0, i)->x);
	}
	fprintf_s(plik, "\n");

	for (int i = 0; i < paramnet->getxsize(); i++)
	{
		fprintf_s(plik, "%.16lf,", paramnet->getnode(i, 0)->t);
		for (int j = 0; j < paramnet->getysize(); j++)
		{
			fprintf_s(plik, "%.16lf,", paramnet->getnode(i, j)->value);
		}
		fprintf_s(plik, "\n");
	}

	fclose(plik);
}

net* analityczna(double h, double deltat)
{

	int hcount = (2. * a) / h;
	int deltatcount = tmax / deltat;

	net* localnet = new net(hcount, h, deltatcount, deltat);

	//rozwiązanie
	for (int i = 0; i < localnet->getxsize(); i++)
	{
		for (int j = 0; j < localnet->getysize(); j++)
		{
			auto tnode = localnet->getnode(i, j);
			tnode->value = rozw_analityczne(tnode->x, tnode->t);
		}

	}

	return localnet;

}

//TODO: Poprawić
net* crankanicolson_thomasa(double h, double deltat)
{
	const double lambda = D * (deltat / (h * h));
	assert(lambda < 1.2);
	assert(lambda > 0.8);


	int hcount = (2. * a) / h;
	int deltatcount = tmax / deltat;

	//inicjalizacja
	net* localnet = new net(hcount, h, deltatcount, deltat);
	for (int i = 0; i < localnet->getysize(); i++)
	{
		auto tnode = localnet->getnode(0, i);
		tnode->value = warunek_poczatkowy(tnode->x);
	}

	vector<double> up;
	vector<double> diag;
	vector<double> low;
	up.reserve(localnet->getysize() - 1);
	diag.reserve(localnet->getysize());
	low.reserve(localnet->getysize() - 1);

	up.push_back(0);
	for (int i = 1; i < localnet->getysize()-1; i++)
	{
		up.push_back(lambda / 2);
		low.push_back(lambda / 2);
	}
	low.push_back(0);

	diag.push_back(1);
	for (int i = 1; i < localnet->getysize()-1; i++)
	{
		diag.push_back(-(1 + lambda));
	}
	diag.push_back(1);
	
	Tridiagonal tri(low, diag, up);
	DU du(tri.siz);
	tri.ThomasMatrix(du);


	//rozwiązanie
	for (int i = 1; i < localnet->getxsize(); i++)
	{

		vector<double> b;
		b.reserve(localnet->getysize());
		//b.push_back(0);
		{
			auto node2 = localnet->getnode(i - 1, 0);
			auto node3 = localnet->getnode(i - 1, 1);
			b.push_back(-(((1 - lambda) * node2->value) + ((lambda / 2) * node3->value)));
		}
		for (int j = 1; j < localnet->getysize() - 1; j++)
		{
			auto node1 = localnet->getnode(i - 1, j - 1);
			auto node2 = localnet->getnode(i - 1, j);
			auto node3 = localnet->getnode(i - 1, j + 1);

			b.push_back(-(((lambda / 2) * node1->value) + ((1 - lambda) * node2->value) + ((lambda / 2) * node3->value)));
		}

		//b.push_back(0);
		{
			auto node1 = localnet->getnode(i - 1, localnet->getysize() - 2);
			auto node2 = localnet->getnode(i - 1, localnet->getysize() - 1);
			b.push_back(-(((lambda / 2) * node1->value) + ((1 - lambda) * node2->value)));
		}
		node_solver(tri, du, b, localnet->body[i]);

	}

	return localnet;
}

net* crankanicolson_gaussa_seidela(double h, double deltat)
{
	const double lambda = D * (deltat / (h * h));
	assert(lambda < 1.2);
	assert(lambda > 0.8);


	int hcount = (2. * a) / h;
	int deltatcount = tmax / deltat;

	//inicjalizacja
	net* localnet = new net(hcount, h, deltatcount, deltat);
	for (int i = 0; i < localnet->getysize(); i++)
	{
		auto tnode = localnet->getnode(0, i);
		tnode->value = warunek_poczatkowy(tnode->x);
	}

	//TUpper U(localnet->getysize());
	FullMatrix LU(localnet->getysize());
	LU.zapisz(0, 1) = 0;
	for (int i = 1; i < localnet->getysize() - 1; i++)
	{
		LU.zapisz(i, i + 1) = lambda / 2;
	}
	for (int i = 1; i < localnet->getysize() - 1; i++)
	{
		LU.zapisz(i, i - 1) = lambda / 2;
	}
	LU.zapisz(localnet->getysize()-1, localnet->getysize() - 2) = 0;
	Diagonal dinv(localnet->getysize());
	for (int i = 1; i < localnet->getysize()-1; i++)
	{
		dinv.zapisz(i, i) = -(1 + lambda);
	}
	dinv.zapisz(0, 0) = 1;
	dinv.zapisz(localnet->getysize()-1, localnet->getysize()-1) = 1;
	dinv.Invert();



	//rozwiązanie
	for (int i = 1; i < localnet->getxsize(); i++)
	{
		Wector x(localnet->getysize());
		for ( int j = 0; j < localnet->getysize(); j++)
		{
			x.zapisz(j) = localnet->getnode(i - 1, j)->value;
		}
		unsigned int n = -1;
		Wector x1(x);
		Wector E(localnet->getysize());
		Wector F(localnet->getysize());

		unsigned int iter = 400;
		double ftol = 1e-15;
		double etol = 1e-15;

		Wector b(localnet->getysize());
		{
			auto node2 = localnet->getnode(i - 1, 0);
			auto node3 = localnet->getnode(i - 1, 1);
			b.zapisz(0) = -(((1 - lambda) * node2->value) + ((lambda / 2) * node3->value));
		}
		for (int j = 1; j < localnet->getysize() - 1; j++)
		{
			auto node1 = localnet->getnode(i - 1, j - 1);
			auto node2 = localnet->getnode(i - 1, j);
			auto node3 = localnet->getnode(i - 1, j + 1);

			b.zapisz(j) = -(((lambda / 2) * node1->value) + ((1 - lambda) * node2->value) + ((lambda / 2) * node3->value));
		}
		{
			auto node1 = localnet->getnode(i - 1, localnet->getysize() - 2);
			auto node2 = localnet->getnode(i - 1, localnet->getysize() - 1);
			b.zapisz(localnet->getysize() -1) = -(((lambda / 2) * node1->value) + ((1 - lambda) * node2->value));
		}

		do
		{
			for (unsigned int k = 0; k < x.siz; k++)
			{
				double temp = 0;
				for (unsigned int l = 0; l < x.siz; l++)
					temp += LU.at(k, l) * x.at(l);
				x.zapisz(k) = dinv.at(k, k) * (b.at(k) - temp);
			}
			x.SubAbs(x1, E);
			LU.Multiply(x, F);
			F.SubAbs(b, F);
			n++;
			x1 = x;
		} while ((n < iter) && (E.Max() > etol) && (F.Max() > ftol));
		
		for (int j = 0; j < localnet->getysize(); j++)
		{
			localnet->getnode(i, j)->value = x.at(j);
		}
		if (i % 20 == 0)
		{
			printf_s("Ukonczono w %lf%%\n", (((double) i / localnet->getxsize()) * 100));
		}
	}

	return localnet;
}

int main()
{
	double h = 0.05;
	//double deltat = h * h * 0.4;
	double deltat = h * h;

	cout << "start" << endl;
	/*auto bezp = bezposrednia_eulera(h, deltat);
	dumpnet(bezp, "tempbezp.csv");
	delete bezp;*/

	/*auto cnik_thom = crankanicolson_thomasa(h, deltat);
	dumpnet(cnik_thom, "tempcnikthom.csv");
	delete cnik_thom;*/
	
	auto cnik_seid = crankanicolson_gaussa_seidela(h, deltat);
	dumpnet(cnik_seid, "tempcnikseid.csv");
	delete cnik_seid;

	cout << "bezp" << endl;
	auto anal = analityczna(h, deltat);
	cout << "anal" << endl;

	
	dumpnet(anal, "tempanal.csv");

	
	delete anal;

	return 0;
}
