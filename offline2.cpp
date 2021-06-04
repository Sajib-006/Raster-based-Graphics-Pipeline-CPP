#include<bits/stdc++.H>
#include<bits/stdc++.H>
#include "bitmap_image.hpp"

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
struct Triangle{
    Point3D points[3];
    int color[3];
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
void printTriangle(struct Triangle t)
{
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++) printf("%f ",t.points[i].arr[j]);
        printf("\n");
    }
    for(int i=0; i<3; i++) printf("%d ",t.color[i]);
    printf("\n\n");

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
double make_line_eqn(Point3D p1, Point3D p2, double y)
{
    double x,x1,x2,y1,y2;
    x1 = p1.arr[0];
    y1 = p1.arr[1];
    x2 = p2.arr[0];
    y2 = p2.arr[1];
    x = ((x1 - x2)/(y1 - y2)*(y - y1)) + x1;
    return x;
}
double calculate_z_value(struct Triangle t, double ys, double xa, double xb, double xp)
{
    double z1,z2,z3,y1,y2,y3,za,zb,zp;
    z1 = t.points[0].arr[2];
    z2 = t.points[1].arr[2];
    z3 = t.points[2].arr[2];
    y1 = t.points[0].arr[1];
    y2 = t.points[1].arr[1];
    y3 = t.points[2].arr[1];
    za = z1 - ((z1 - z2)*(y1 - ys)/(y1 - y2));
    zb = z1 - ((z1 - z3)*(y1 - ys)/(y1 - y3));
    zp = zb - ((zb - za)*(xb - xp)/(xb - xa));
    return zp;

}
void scan_conversion(ifstream& fin, ifstream& config, ofstream& fout)
{
    double screen_width, screen_height;
    double neg_x_limit, pos_x_limit, bottom_y_limit, top_y_limit, front_z_limit, rear_z_limit;
    // read config file
    config >> screen_width >> screen_height;
    config >> neg_x_limit;
    config >> bottom_y_limit;
    config >> front_z_limit >> rear_z_limit;

    pos_x_limit = neg_x_limit * -1;
    top_y_limit = bottom_y_limit * -1;
    cout<<"Stage4 output"<<endl<<endl;
    cout<< front_z_limit<<" "<< rear_z_limit<<" "<<neg_x_limit<<" "<<pos_x_limit<<" "<<bottom_y_limit<<" "<<top_y_limit<<endl;
    cout<<triangle_cnt<<endl;

    struct Triangle triangles[triangle_cnt];
    int cnt = triangle_cnt;
    int idx = 0;
    int min_v = 0, max_v = 255;
    while(cnt--)
    {
       for(int i=0; i<3; i++)
       {
           for(int j=0; j<3; j++)
           {
               fin>> triangles[idx].points[i].arr[j];
           }
           triangles[idx].points[i].arr[3] = 1;
       }
       for(int i=0; i<3; i++) triangles[idx].color[i] = rand()%(max_v-min_v + 1) + min_v;
       idx++;
    }
    for(int i=0; i<triangle_cnt; i++) printTriangle(triangles[i]);

    double dx, dy, top_y, left_x, mid_x, mid_y;
    dx = (pos_x_limit - neg_x_limit)/screen_width;
    dy = (top_y_limit - bottom_y_limit)/screen_height;
    top_y = top_y_limit - dy/2;
    left_x = neg_x_limit + dx/2;

    //initialize dynamic z_buffer array
    double** z_buffer = new double*[int(screen_width)];
    for(int i = 0; i < screen_width; ++i)
        z_buffer[i] = new double[int(screen_height)];
    //double z_buffer = new double[screen_width][screen_height];

    //initialize bitmap image
    bitmap_image image(screen_width,screen_height);
    double z_max = rear_z_limit;

    for(int row_no=0; row_no<screen_width; row_no++)
    {
        for(int col_no=0; col_no<screen_height; col_no++)
        {
            z_buffer[row_no][col_no] = z_max;
            image.set_pixel(row_no,col_no,0,0,0);
        }
    }
    image.save_image("4_out/background.bmp");;
    cout<<"after image background set"<<endl;

    //applying procedure
    double z_value,ys,xp,max_x, min_x, top_scanline, bottom_scanline, line1_x, line2_x, line3_x, left_int_col, right_int_col;
    int top_cell, bottom_cell,left_cell,right_cell;
    struct Triangle t;

    for(int i=0; i<triangle_cnt; i++)
    {
        cout<<"Triangle "<<i<<endl;
        t = triangles[i];
        top_scanline = max(max(t.points[0].arr[1],t.points[1].arr[1]),t.points[2].arr[1]);
        bottom_scanline = min(min(t.points[0].arr[1],t.points[1].arr[1]),t.points[2].arr[1]);
        //top_scanline -= dy/2;
        if(top_scanline > top_y) top_scanline = top_y;
        if(bottom_scanline < bottom_y_limit) bottom_scanline = bottom_y_limit;
        top_cell = (top_y - top_scanline) / dy;
        bottom_cell = (top_y - bottom_scanline) / dy;

        for(int row_no=top_cell; row_no<bottom_cell; row_no++)
        {
            cout<<"row "<<row_no<<endl;
            ys = row_no - dy/2;
            line1_x = make_line_eqn(t.points[0],t.points[1],ys);
            line2_x = make_line_eqn(t.points[1],t.points[2],ys);
            line3_x = make_line_eqn(t.points[0],t.points[2],ys);

            if(line1_x < left_x || line1_x > pos_x_limit)
            {
                left_int_col = min(line2_x, line3_x);
                right_int_col = max(line2_x, line3_x);
            }
            else if(line2_x < left_x || line2_x > pos_x_limit)
            {
                left_int_col = min(line1_x, line3_x);
                right_int_col = max(line1_x, line3_x);
            }
            else if(line3_x < left_x || line3_x > pos_x_limit)
            {
                left_int_col = min(line1_x, line2_x);
                right_int_col = max(line1_x, line2_x);
            }
            //left_int_col = min(min(line1_x,line2_x),line3_x);
            //right_int_col = max(max(line1_x,line2_x),line3_x);
            //left_int_col += dx/2;

            //if(left_int_col < left_x) left_int_col = left_x;
            //if(right_int_col > pos_x_limit) right_int_col = pos_x_limit;

            left_cell = (left_int_col - left_x) / dx;
            right_cell = (right_int_col - left_x) / dx;
            for(int col_no=left_cell; col_no<right_cell; col_no++)
            {
                //mid_y = top_y - row_no * dy;
                //mid_x = left_x + col_no * dx;
                xp = col_no + dx/2;
                z_value = calculate_z_value(t, ys , left_int_col, right_int_col, xp);
                if(z_value >= front_z_limit && z_value < z_buffer[row_no][col_no]) {
                    z_buffer[row_no][col_no] = z_value;
                    image.set_pixel(row_no,col_no,t.color[0],t.color[1],t.color[2]);
                    cout<<"Inside "<<i<<" "<<row_no<<" "<<col_no<<endl;
                }


            }//end col-loop
        }//end row-loop
    }//end triangle-loop
    image.save_image("4_out/output.bmp");;
    cout<<"done"<<endl;
    for(int i=0; i<screen_width; i++){
        for(int j=0; j<screen_height; j++){
            if(z_buffer[i][j]< z_max)
                fout<<z_buffer[i][j]<<" ";
        }
        fout<<endl;
    }


    //free z_buffer memory
    for(int i = 0; i < screen_width; ++i) {
        delete [] z_buffer[i];
    }
    //free image memory



}
int main()
{
    ifstream fin1( "4/scene.txt" );
    ofstream fout1("4_out/stage1.txt");
    ifstream fin2("4_out/stage1.txt");
    ofstream fout2("4_out/stage2.txt");
    ifstream fin3("4_out/stage2.txt");
    ofstream fout3("4_out/stage3.txt");
    ifstream fin4("4_out/stage3.txt");
    ofstream fout4("4_out/z_buffer.txt");
    ifstream config("4/config.txt");

    modeling_transformation(fin1,fout1);
    view_transformation(fin2,fout2);
    project_transformation(fin3,fout3);
    scan_conversion(fin4,config,fout4);
    /*double screen_width =500, screen_height =500;
    bitmap_image image(screen_width,screen_height);

    for(int row_no=0; row_no<screen_width; row_no++)
    {
        for(int col_no=0; col_no<screen_height; col_no++)
        {

            image.set_pixel(row_no,col_no,0,0,0);
        }
    }
    image.save_image("bc.bmp");;
    cout<<"after image background set"<<endl;*/
}
