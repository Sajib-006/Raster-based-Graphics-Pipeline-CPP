#include<bits/stdc++.H>
using namespace std;

struct Point
{
	double x = 0.0,y = 0.0,z = 0.0, w = 1.0;
};
struct Vector
{
    double x = 0.0,y = 0.0,z = 0.0, w = 1.0;
};
struct Matrix
{
    /*double x1 = 1.0, x2 = 0.0, x3 = 0.0, x4 = 0.0;
    double y1 = 0.0, y2 = 1.0, y3 = 0.0, y4 = 0.0;
    double z1 = 0.0, z2 = 0.0, z3 = 1.0, z4 = 0.0;
    double w1 = 0.0, w2 = 0.0, w3 = 0.0, w4 = 1.0;*/
    double arr[4][4];
};
stack<Matrix> s;
struct Point eye;
struct Point look;
struct Point up;
double fovY, aspectRatio, near, far;

void printPoint(struct Point p)
{
    printf("%f %f %f %f\n",p.x,p.y,p.z,p.w);
}
void printMatrix(struct Matrix m)
{
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++) printf("%f ",m.arr[i][j]);
        printf("\n");
    }
    printf("\n");

}
struct Matrix initializeMatrix()
{
    struct Matrix m;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++) m.arr[i][j] = 0.0;
    m.arr[3][3] = 1.0;
    return m;
}
struct Matrix matrixMultiplication(struct Matrix matrix1, struct Matrix matrix2)
{
    struct Matrix result;
    for(int i=0; i<4; i++){
       for(int j=0; j<4; j++){
            result.arr[i][j] = 0.0;
            for(int k=0; k<4; k++)
                result.arr[i][j] += matrix1.arr[i][k] * matrix2.arr[k][j];
        }
    }
    return result;

}
struct Matrix CrossMultiplication(struct Matrix m, struct Vector v)
{



}

struct Vector normalize(struct Vector a)
{
    double val = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    a.x /= val;
    a.y /= val;
    a.z /= val;
    return a;
}
struct Matrix genIdentityMatrix()
{
    struct Matrix m;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            if(i==j) m.arr[i][j] = 1.0;
            else m.arr[i][j] = 0.0;
    return m;
}
struct Matrix genTranslationMatrix(double x, double y, double z)
{
    struct Matrix m;
    m = initializeMatrix();
    m.arr[0][3] = x;
    m.arr[1][3] = y;
    m.arr[2][3] = z;
    return m;
}
struct Matrix genScaleMatrix(double x, double y, double z)
{
    struct Matrix m;
    m = initializeMatrix();
    m.arr[0][0] = x;
    m.arr[1][1] = y;
    m.arr[2][2] = z;
    return m;
}
struct Matrix genRotationMatrixUtil(struct Vector x,double angle, struct Vector a)
{

}
struct Matrix genRotationMatrix(double angle, struct Vector a)
{

}

int main()
{
    cout<<"Hello"<<endl;
    struct Matrix identity_matrix = genIdentityMatrix();
    printMatrix(identity_matrix);
    s.push(identity_matrix);
    ifstream fin( "1/scene.txt" );

    fin >> eye.x >> eye.y >> eye.z;
    fin >> look.x >> look.y >> look.z;
    fin >> up.x >> up.y >> up.z;
    fin >> fovY >> aspectRatio >> near >> far;

    printPoint(eye);
    //while(getline( fin, line ))
    string command;
    while(true)
    {
        fin >> command;
        if(command == "triangle"){
           cout<<"T"<<endl;
           struct Matrix triangle = initializeMatrix();
           for(int i=0; i<3; i++){
              for(int j=0; j<3; j++){
                fin>> triangle.arr[i][j];
              }
           }

           printMatrix(triangle);
           printMatrix(s.top());
           //printMatrix(matrixMultiplication(triangle,s.top()));
           s.push(matrixMultiplication(triangle,s.top()));
           printMatrix(s.top());

        }
        else if(command == "translate"){
            cout<<"Trans"<<endl;
            double x,y,z;
            fin>> x >> y >> z;
            struct Matrix T = genTranslationMatrix(x,y,z);
            s.push(matrixMultiplication(T,s.top()));
        }
        else if(command == "scale"){
            cout<<"S"<<endl;
            double x,y,z;
            fin>> x >> y >> z;
            struct Matrix T = genScaleMatrix(x,y,z);
            s.push(matrixMultiplication(T,s.top()));
        }
        else if(command == "rotate"){
            cout<<"R"<<endl;
            struct Vector a;
            double angle;
            fin>> angle >> a.x >> a.y >> a.z;
            a = normalize(a);
        }
        else if(command == "push"){
            cout<<"Push"<<endl;
        }
        else if(command == "pop"){
            cout<<"pop"<<endl;
        }
        else if(command == "end") break;
    }

}
