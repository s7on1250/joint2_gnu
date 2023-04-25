#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdio>
#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif
using namespace std;
class Matrix {
protected:
    double **matrix;
    int rows = 0;
    int columns = 0;
    bool isCorrupted = false;
public:
    Matrix() {}

    Matrix(int n, int m) {
        if (rows != 0 || columns != 0) {
            for (int i = 0; i < rows; i++) {
                delete matrix[i];
            }
            delete matrix;
        }
        this->rows = n;
        this->columns = m;

        matrix = new double *[n];
        for (int i = 0; i < n; i++) {
            matrix[i] = new double[m];

        }
    }

    Matrix(const Matrix &another) : Matrix(another.rows, another.columns) {
        this->isCorrupted = another.isCorrupted;
        if (another.rows >= 0 && another.columns >= 0) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    this->matrix[i][j] = another.matrix[i][j];

                }
            }
        }

    }

    double getElement(int i, int j) {
        if (i >= 0 && i < getRows() && j >= 0 && j < getColumns()) {
            double element = matrix[i][j];
            return element;
        }
        return 0;
    }

    void setElement(int i, int j, double num) {
        if (i >= 0 && i < getRows() && j >= 0 && j < getColumns()) {
            matrix[i][j] = num;
        }
    }

    int getRows() const {
        return rows;
    }

    int getColumns() const {
        return columns;
    }

    bool getCorruption() const {
        return isCorrupted;
    }

    friend istream &operator>>(istream &in, Matrix &m) {
        int R, C;
        in >> R;
        in >> C;
        m = *new Matrix(R, C);
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < C; j++) {
                in >> m.matrix[i][j];
            }
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, Matrix &m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                double element = m.matrix[i][j];
                out << element << " ";
            }
            out << endl;
        }
        return out;
    }

    Matrix &operator=(const Matrix &another) = default;

    Matrix &operator+(const Matrix &another) {
        auto *temp = new Matrix(*this);
        if (rows != another.rows || columns != another.columns) {
            cout << "Error: the dimensional problem occurred" << endl;
            temp->isCorrupted = true;
            return *temp;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp->matrix[i][j] += another.matrix[i][j];
            }
        }
        return *temp;
    }

    Matrix &operator-(const Matrix &another) {
        auto *temp = new Matrix(*this);
        if (rows != another.rows || columns != another.columns) {
            cout << "Error: the dimensional problem occurred" << endl;
            temp->isCorrupted = true;
            return *temp;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp->matrix[i][j] -= another.matrix[i][j];

            }
        }
        return *temp;
    }

    Matrix &operator*(const Matrix &another) {
        auto *tempo = new Matrix(*this);
        if (columns != another.rows) {
            cout << "Error: the dimensional problem occurred" << endl;
            tempo->isCorrupted = true;
            return *tempo;
        }
        Matrix *temp = new Matrix(rows, another.columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < another.columns; j++) {
                temp->matrix[i][j] = 0;
                for (int k = 0; k < another.rows; k++) {
                    temp->matrix[i][j] += matrix[i][k] * another.matrix[k][j];
                }
            }
        }
        return *temp;
    }
    Matrix &operator*(const double &alpha) {

        Matrix *temp = new Matrix(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp->matrix[i][j] = matrix[i][j]*alpha;
            }
        }
        return *temp;
    }
    Matrix transpose() {
        Matrix temp(columns, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp.matrix[j][i] = this->matrix[i][j];
            }
        }
        return temp;
    }

    ~Matrix() {
        for (int i = 0; i < rows; i++) {
            delete matrix[i];
        }
        delete matrix;
    }
};

class SquareMatrix : public Matrix {
protected:
    int size;
public:
    SquareMatrix() = default;

    explicit SquareMatrix(int N) : Matrix(N, N) {
        size = rows;
    }

    SquareMatrix(Matrix m) : SquareMatrix(m.getColumns()) {
        if (m.getColumns() == m.getRows()) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    matrix[i][j] = m.getElement(i, j);
                }
            }
        }
    }

    int getSize() {
        return size;
    }

    friend istream &operator>>(istream &in, SquareMatrix &m) {
        int N;
        in >> N;
        m = *new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                in >> m.matrix[i][j];
            }
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, SquareMatrix &m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                double element = m.matrix[i][j];
                out << element << " ";
            }
            out << endl;
        }
        return out;
    }

    SquareMatrix &operator=(const Matrix &another) {
        new(this)SquareMatrix(another);
        return *this;
    };
    Matrix inverse();
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    this->matrix[i][j] = 1;
                } else
                    this->matrix[i][j] = 0;
            }

        }
    }

    friend istream &operator>>(istream &in, IdentityMatrix &m) {
        int N;
        in >> N;
        m = *new IdentityMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                in >> m.matrix[i][j];
            }
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, IdentityMatrix &m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                double element = m.matrix[i][j];
                out << element << " ";
            }
            out << endl;
        }
        return out;
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(int i, int j, SquareMatrix m) : IdentityMatrix(m.getSize()) {
        i -= 1;
        j -= 1;
        double value = -m.getElement(i, j) / m.getElement(j, j);
        this->matrix[i][j] = value;

    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int i, int j, SquareMatrix m) : IdentityMatrix(m.getSize()) {
        i -= 1;
        j -= 1;
        int temp;
        for (int n = 0; n < m.getSize(); n++) {
            temp = matrix[i][n];
            matrix[i][n] = matrix[j][n];
            matrix[j][n] = temp;
        }
    }
};

Matrix SquareMatrix::inverse() {
    SquareMatrix a = *this;
    SquareMatrix result = IdentityMatrix(size);
    int currentColumn = 0;
    int pivot = 0;
    int step = 0;
    //making U matrix
    while (currentColumn != a.getSize()) {
        double max = INT64_MIN;
        int idx = 0;
        //permutation
        for (int i = pivot + 1; i < a.getSize(); i++) {
            if (abs(a.getElement(i, currentColumn)) > max) {
                max = abs(a.getElement(i, currentColumn));
                idx = i;
            }
        }
        if (max > abs(a.getElement(pivot, currentColumn))) {
            PermutationMatrix permutationMatrix(pivot + 1, idx + 1, a);
            a = permutationMatrix * a;
            result = permutationMatrix * result;

        }
        //elimination (making U matrix)
        for (int i = pivot + 1; i < a.getSize(); i++) {
            if (a.getElement(i, currentColumn) != 0) {
                EliminationMatrix eliminationMatrix(i + 1, currentColumn + 1, a);
                a = eliminationMatrix * a;
                result = eliminationMatrix * result;
            }
        }
        currentColumn++;
        pivot++;
    }
    //making D matrix
    pivot--;
    currentColumn--;
    while (currentColumn != -1) {
        //elimination (making D matrix)
        for (int i = pivot - 1; i >= 0; i--) {
            if (a.getElement(i, currentColumn) != 0) {
                EliminationMatrix eliminationMatrix(i + 1, currentColumn + 1, a);
                a = eliminationMatrix * a;
                result = eliminationMatrix * result;
            }
        }
        currentColumn--;
        pivot--;
    }
    for (int i = 0; i < a.getSize(); i++) {
        double coef = a.getElement(i, i);
        for (int j = 0; j < result.getColumns(); j++) {
            a.setElement(i, j, a.getElement(i, j) / coef);
            result.setElement(i, j, result.getElement(i, j) / coef);
        }
    }
    return result;
}

class ColumnVector : public Matrix{
public:
    ColumnVector() = default;
    ColumnVector(int n) : Matrix(n,1){}
    friend istream &operator>>(istream &in, ColumnVector &m) {
        int N;
        in >> N;
        m = *new ColumnVector(N);
        for (int i = 0; i < N; i++) {
            in >> m.matrix[i][0];
        }
        return in;
    }
    double norm(){
        double sum = 0;
        for(int i=0;i<rows;i++){
            sum+=matrix[i][0];
        }
        return sum/ rows;
    }
    ColumnVector &operator*(const double &alpha) {

        ColumnVector *temp = new ColumnVector(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp->matrix[i][j] = matrix[i][j]*alpha;
            }
        }
        return *temp;
    }

};
int main() {
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME,"w");
#else
    FILE* pipe = popen(GNUPLOT_NAME,"w");
#endif
    int m,n;
    cin>>m;
    double ts[m];
    ColumnVector bv(m);
    for(int i=0;i<m;i++){
        double b;
        cin>>ts[i]>>b;
        bv.setElement(i,0,b);
    }
    cin>>n;
    Matrix a(m,n+1);
    for(int i=0;i<m;i++){
        for(int j=0;j<n+1;j++){
            a.setElement(i,j,pow(ts[i],j));
        }
    }
    cout<<setprecision(4)<<fixed<<"A:"<<endl<<a;
    Matrix a_T = a.transpose();
    SquareMatrix a_T_a = a_T*a;
    cout<<setprecision(4)<<fixed<<"A_T*A:"<<endl<<a_T_a;
    SquareMatrix a_T_a_I = a_T_a.inverse();
    cout<<setprecision(4)<<fixed<<"(A_T*A)^-1:"<<endl<<a_T_a_I;
    Matrix a_T_b = a_T*bv;
    cout<<setprecision(4)<<fixed<<"A_T*b:"<<endl<<a_T_b;
    Matrix result = a_T_a_I*a_T_b;
    cout<<setprecision(4)<<fixed<<"x~:"<<endl<<result;

    //gnu
    //fprintf(pipe,"%s\n","set term wxt");
    fprintf(pipe,"%s\n","set multiplot");
    fprintf(pipe,"%s\n","set xrange [-10:10]");
    fprintf(pipe,"%s\n","set yrange [-10:10]");
    //initial points
    fprintf(pipe,"%s\n","plot '-' using 1:2 title 'initial points' with points");
    for(int i=0;i<m;i++){
        fprintf(pipe,"%f\t%f\n",ts[i],bv.getElement(i,0));
    }
    fprintf(pipe,"%s\n","e");
    //approximation graph
    fprintf(pipe,"%s","plot ");
    for(int i=0;i<=n;i++){
        fprintf(pipe,"(%f)*x**%d",result.getElement(n-i,0),n-i);
        if(i!=n){
            fprintf(pipe,"+");
        }
    }
    fprintf(pipe,"%s\n","using 1:2 title 'approximation graph' with lines ");
    fprintf(pipe,"%s\n","unset multiplot");
    fflush(pipe);
    return 0;
}
