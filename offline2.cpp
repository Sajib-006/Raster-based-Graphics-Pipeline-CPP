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
double pi = 3.1416;
int triangle_cnt = 0;
stack<Matrix> s;
stack<string> command_stack;
/*struct Point eye;
struct Point look;
struct Point up;*/
struct Vector eye,look,up;
double fovY, aspectRatio, near, far;

void printPoint(struct Point p)
{
    printf("%f %f %f %f\n",p.x,p.y,p.z,p.w);
}
void printPoint3D(struct Point3D p,ofstream& fout)
{
    if(p.arr[3]!=1){
        p.arr[0]/=p.arr[3];
        p.arr[1]/=p.arr[3];
        p.arr[2]/=p.arr[3];
        p.arr[3]/=p.arr[3];
    }
    printf("%f %f %f %f\n",p.arr[0],p.arr[1],p.arr[2],p.arr[3]);
    fout <<p.arr[0]<<" "<<p.arr[1]<<" "<<p.arr[2]<<endl;
}
void printVector(struct Vector v)
{
    printf("%f %f %f %f\n",v.x,v.y,v.z,v.w);
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
struct Vector VecCrossMul(struct Vector a, struct Vector b)
{
    struct Vector res;
    res.x = a.y * b.z - a.z * b.y;
    res.y = a.z * b.x - a.x * b.z;
    res.z = a.x * b.y - a.y * b.x;
    res.w = 1;
    return res;
}
struct Vector VecSub(struct Vector a, struct Vector b)
{
    struct Vector res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    res.w = 1;
    return res;
}
struct Vector VecAdd(struct Vector a, struct Vector b)
{
    struct Vector res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    res.w = 1;
    return res;
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
    v3.x = a.y*v.z - a.z*v.y;
    v3.y = a.z*v.x - a.x*v.z;
    v3.z = a.x*v.y - a.y*v.x;
    v3 = ScalarVectorMul(sin(angle),v3);

    if(angle<0.0) angle*=-1;
    v1 = ScalarVectorMul(cos(angle),v);
    double val = (a.x*v.x + a.y*v.y + a.z*v.z)*(1.0-cos(angle));
    v2 = ScalarVectorMul(val,a);

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
    i.x = 1.0;
    j.y = 1.0;
    k.z = 1.0;
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
void modeling_transformation(ifstream& fin, ofstream& fout)
{
    struct Matrix identity_matrix = genIdentityMatrix();
    s.push(identity_matrix);

    fin >> eye.x >> eye.y >> eye.z;
    fin >> look.x >> look.y >> look.z;
    fin >> up.x >> up.y >> up.z;
    fin >> fovY >> aspectRatio >> near >> far;

    cout<<"Stage1 output\n--------------\n";
    string command;
    while(true)
    {
        fin >> command;
        if(command == "triangle"){
           triangle_cnt++;
           Point3D p1,p2,p3;
           for(int i=0; i<3; i++) fin>> p1.arr[i];
           for(int i=0; i<3; i++) fin>> p2.arr[i];
           for(int i=0; i<3; i++) fin>> p3.arr[i];
           p1.arr[3] = p2.arr[3] = p3.arr[3] = 1;

           printPoint3D(matrixPointMul(s.top(),p1),fout);
           printPoint3D(matrixPointMul(s.top(),p2),fout);
           printPoint3D(matrixPointMul(s.top(),p3),fout);
           cout<<endl;
           fout<<endl;

        }
        else if(command == "translate"){
            command_stack.push("T");
            //cout<<"Trans"<<endl;
            double x,y,z;
            fin>> x >> y >> z;
            struct Matrix T = genTranslationMatrix(x,y,z);
            s.push(matrixMultiplication(s.top(),T));
        }
        else if(command == "scale"){
            command_stack.push("S");
            //cout<<"S"<<endl;
            double x,y,z;
            fin>> x >> y >> z;
            struct Matrix T = genScaleMatrix(x,y,z);
            s.push(matrixMultiplication(s.top(),T));
        }
        else if(command == "rotate"){
            command_stack.push("R");
            //cout<<"R"<<endl;
            struct Vector a;
            double angle;

            fin>> angle >> a.x >> a.y >> a.z;
            a = normalize(a);

            //cout<<"angle: "<<angle<<endl;
            struct Matrix R = genRotationMatrix(angle*pi/180.0,a);
            //printMatrix(R);
            s.push(matrixMultiplication(s.top(),R));
        }
        else if(command == "push"){
            command_stack.push("push");
            //cout<<"Push"<<endl;
        }
        else if(command == "pop"){
            //cout<<"pop"<<endl;
            while(command_stack.top()!="push"){
                command_stack.pop();
                s.pop();
            }
            command_stack.pop();
        }
        else if(command == "end") break;
    }
    fin.close();
    fout.close();
}

struct Matrix view_transformation_util()
{
    struct Vector l,r,u;
    struct Matrix V, R, T;

    l = VecSub(look,eye);
    l = normalize(l);
    r = VecCrossMul(l,up);
    r = normalize(r);
    u = VecCrossMul(r,l);

    T = genIdentityMatrix();
    T.arr[0][3] = - eye.x;
    T.arr[1][3] = - eye.y;
    T.arr[2][3] = - eye.z;

    R = initializeMatrix();
    R.arr[0][0] = r.x;
    R.arr[0][1] = r.y;
    R.arr[0][2] = r.z;
    R.arr[1][0] = u.x;
    R.arr[1][1] = u.y;
    R.arr[1][2] = u.z;
    R.arr[2][0] = -l.x;
    R.arr[2][1] = -l.y;
    R.arr[2][2] = -l.z;
    V = matrixMultiplication(R,T);
    return V;
}
void view_transformation(ifstream& fin, ofstream& fout)
{
    int cnt;
    struct Matrix V;
    V = view_transformation_util();
    cout<<"Stage2 output\n--------------\n";
    cnt = triangle_cnt;
    //cout<<triangle_cnt<<endl;
    while(cnt--)
    {
        //cout<<cnt<<endl;
       Point3D p1,p2,p3;
       for(int i=0; i<3; i++) fin>> p1.arr[i];
       for(int i=0; i<3; i++) fin>> p2.arr[i];
       for(int i=0; i<3; i++) fin>> p3.arr[i];
       p1.arr[3] = p2.arr[3] = p3.arr[3] = 1;

       printPoint3D(matrixPointMul(V,p1),fout);
       printPoint3D(matrixPointMul(V,p2),fout);
       printPoint3D(matrixPointMul(V,p3),fout);
       cout<<endl;
       fout<<endl;
    }
    fin.close();
    fout.close();

}
struct Matrix project_transformation_util()
{
    double fovX,t,r;
    struct Matrix P;
    fovX = fovY * aspectRatio;
    t = near * tan(fovY/2.0*pi/180.0);
    r = near * tan(fovX/2.0*pi/180.0);
    P = initializeMatrix();
    P.arr[0][0] = near/r;
    P.arr[1][1] = near/t;
    P.arr[2][2] = -(far+near)/(far-near);
    P.arr[2][3] = -(2*far*near)/(far-near);
    P.arr[3][2] = -1;
    return P;

}
void project_transformation(ifstream& fin, ofstream& fout)
{
    int cnt;
    struct Matrix P;
    P = project_transformation_util();
    cout<<"Stage3 output\n--------------\n";
    cnt = triangle_cnt;
    while(cnt--)
    {
       Point3D p1,p2,p3;
       for(int i=0; i<3; i++) fin>> p1.arr[i];
       for(int i=0; i<3; i++) fin>> p2.arr[i];
       for(int i=0; i<3; i++) fin>> p3.arr[i];
       p1.arr[3] = p2.arr[3] = p3.arr[3] = 1;

       printPoint3D(matrixPointMul(P,p1),fout);
       printPoint3D(matrixPointMul(P,p2),fout);
       printPoint3D(matrixPointMul(P,p3),fout);
       cout<<endl;
       fout<<endl;
    }
    fin.close();
    fout.close();
}
int main()
{
    ifstream fin1( "4/scene.txt" );
    ofstream fout1("4_out/stage1.txt");
    ifstream fin2("4_out/stage1.txt");
    ofstream fout2("4_out/stage2.txt");
    ifstream fin3("4_out/stage2.txt");
    ofstream fout3("4_out/stage3.txt");

    modeling_transformation(fin1,fout1);
    view_transformation(fin2,fout2);
    project_transformation(fin3,fout3);
}
