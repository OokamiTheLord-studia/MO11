#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

using namespace std;

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

int main()
{



}
