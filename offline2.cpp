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
struct Point3D
{
    double arr[4];
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
stack<string> command_stack;
struct Point eye;
struct Point look;
struct Point up;
double fovY, aspectRatio, near, far;

void printPoint(struct Point p)
{
    printf("%f %f %f %f\n",p.x,p.y,p.z,p.w);
}
void printPoint3D(struct Point3D p)
{
    printf("%f %f %f\n",p.arr[0],p.arr[1],p.arr[2]);
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
struct Point3D matrixPointMul(struct Matrix matrix1, struct Point3D p)
{
    //cout<<"In point mul"<<endl;
    struct Point3D result;
    for(int i=0; i<4; i++){
        result.arr[i] = 0.0;
        for(int k=0; k<4; k++)
            result.arr[i] += matrix1.arr[i][k] * p.arr[k];

    }
    //printPoint3D(result);
    return result;

}
struct Matrix matrixMultiplication(struct Matrix matrix1, struct Matrix matrix2)
{
    //cout<<"In mat mul"<<endl;
    struct Matrix result;
    for(int i=0; i<4; i++){
       for(int j=0; j<4; j++){
            result.arr[i][j] = 0.0;
            for(int k=0; k<4; k++)
                result.arr[i][j] += matrix1.arr[i][k] * matrix2.arr[k][j];
        }
    }
    //printMatrix(result);
    return result;

}
struct Vector ScalarVectorMul(double scalar, struct Vector v)
{
    struct Vector v1;
    v1.x = scalar * v.x;
    v1.y = scalar * v.y;
    v1.z = scalar * v.z;
    return v1;
}
struct Matrix ScalarMatrixMul(double scalar, struct Matrix m)
{
    struct Matrix m1;
    //for(int i=0; i<)
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

    m.arr[0][0] = 1;
    m.arr[1][1] = 1;
    m.arr[2][2] = 1;
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
struct Vector RotationUtil(struct Vector v, struct Vector a, double angle)
{
    struct Vector v1,v2,v3,res;
    v1 = ScalarVectorMul(cos(angle),v);
    double val = (a.x*v.x + a.y*v.y + a.z*v.z)*(1-cos(angle));
    v2 = ScalarVectorMul(val,a);
    v3.x = a.y*v.z - a.z*v.y;
    v3.y = a.z*v.x - a.x*v.z;
    v3.z = a.x*v.y - a.y*v.x;
    res.x = v1.x + v2.x + v3.x;
    res.y = v1.y + v2.y + v3.y;
    res.z = v1.z + v2.z + v3.z;
    if(res.x>-0.001 && res.x<0.0) res.x = 0.0;
    if(res.y>-0.001 && res.y<0.0) res.y = 0.0;
    if(res.z>-0.001 && res.z<0.0) res.z = 0.0;
    return res;
}
struct Matrix genRotationMatrix(double angle, struct Vector a)
{
    struct Vector i,j,k,c1,c2,c3;
    struct Matrix res;
    i.x = 1;
    j.y = 1;
    k.z = 1;
    c1 = RotationUtil(i,a,angle);
    c2 = RotationUtil(j,a,angle);
    c3 = RotationUtil(k,a,angle);
    res.arr[0][0] = c1.x;
    res.arr[0][1] = c2.x;
    res.arr[0][2] = c3.x;
    res.arr[1][0] = c1.y;
    res.arr[1][1] = c2.y;
    res.arr[1][2] = c3.y;
    res.arr[2][0] = c1.z;
    res.arr[2][1] = c2.z;
    res.arr[2][2] = c3.z;
    res.arr[3][3] = 1;
    res.arr[0][3] = res.arr[1][3] = res.arr[2][3] = res.arr[3][0] = res.arr[3][1] = res.arr[3][2] = 0;
    return res;
}

int main()
{
    cout<<"Hello"<<endl;
    struct Matrix identity_matrix = genIdentityMatrix();
    printMatrix(identity_matrix);
    s.push(identity_matrix);
    ifstream fin( "3/scene.txt" );

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
           cout<<"Triangle"<<endl;
           struct Matrix triangle = initializeMatrix();
           Point3D p1,p2,p3;
           for(int i=0; i<3; i++) fin>> p1.arr[i];
           for(int i=0; i<3; i++) fin>> p2.arr[i];
           for(int i=0; i<3; i++) fin>> p3.arr[i];
           p1.arr[3] = p2.arr[3] = p3.arr[3] = 1;

           //printMatrix(s.top());
           //for(int i=0;i<command_stack.size();i++) cout<<command_stack[i]<<" ";
           printf("\nPrinting Transformed Point\n");
           printPoint3D(matrixPointMul(s.top(),p1));
           printPoint3D(matrixPointMul(s.top(),p2));
           printPoint3D(matrixPointMul(s.top(),p3));

        }
        else if(command == "translate"){
            command_stack.push("T");
            cout<<"Trans"<<endl;
            double x,y,z;
            fin>> x >> y >> z;
            struct Matrix T = genTranslationMatrix(x,y,z);
            s.push(matrixMultiplication(s.top(),T));
        }
        else if(command == "scale"){
            command_stack.push("S");
            cout<<"S"<<endl;
            double x,y,z;
            fin>> x >> y >> z;
            struct Matrix T = genScaleMatrix(x,y,z);
            s.push(matrixMultiplication(s.top(),T));
        }
        else if(command == "rotate"){
            command_stack.push("R");
            cout<<"R"<<endl;
            struct Vector a;
            double angle;
            double pi = 22.0/7;
            fin>> angle >> a.x >> a.y >> a.z;
            a = normalize(a);
            if(angle<0.0) angle*=-1;
            cout<<"angle: "<<angle<<endl;
            struct Matrix R = genRotationMatrix(angle*pi/180.0,a);
            printMatrix(R);
            s.push(matrixMultiplication(s.top(),R));
        }
        else if(command == "push"){
            command_stack.push("push");
            cout<<"Push"<<endl;
        }
        else if(command == "pop"){
            cout<<"pop"<<endl;
            while(command_stack.top()!="push"){
                command_stack.pop();
                s.pop();
            }
            command_stack.pop();
        }
        else if(command == "end") break;
    }

}
