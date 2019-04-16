#ifndef UTIL_H
#define UTIL_H
#include <QVector>
#include <QDebug>
#include <QVector3D>
#include <QImage>
#include <cmath>
class Point{
public:
    double x,y,z;
    Point(){}
    Point(double x1,double y1,double z1){x = x1;y = y1;z = z1;}
    void init(double x1,double y1,double z1){x = x1;y = y1;z = z1;}
    bool operator ==(const Point &p)
    {
        return fabs(x-p.x)<= 1e-6&&fabs(y-p.y)<= 1e-6&&fabs(z-p.z)<= 1e-6;
    }
    Point operator +(const Point& p1)
    {
        return Point(p1.x+x,p1.y+y,p1.z+z);
    }
    Point operator -(const Point& p1)
    {
        return Point(x-p1.x,y-p1.y,z-p1.z);
    }

    void print()
    {
        qDebug()<<"("<<x<<","<<y<<","<<z<<")";
    }


};
class vec3{  //向量
public:
    double x,y,z;
    vec3(){}
    vec3(double x1,double y1,double z1=0){x = x1;y = y1;z = z1;}
    vec3(Point a){ x = a.x; y = a.y;z = a.z;}
    vec3(Point a,Point b){ x = b.x-a.x; y = b.y-a.y; z = b.z-a.z; }  //a指向b的向量
    void init(double x1,double y1,double z1=0){x = x1;y = y1;z = z1;}
    bool operator ==(const vec3 &p)
    {
        return fabs(x-p.x)<= 1e-6&&fabs(y-p.y)<= 1e-6&&fabs(z-p.z)<= 1e-6;
    }
    vec3 operator +(const vec3& p1)
    {
        return vec3(p1.x+x,p1.y+y,p1.z+z);
    }
    vec3 operator -(const vec3& p1)
    {
        return vec3(x-p1.x,y-p1.y,z-p1.z);
    }
    vec3 operator *(double num)
    {
        return vec3(num*x,num*y,num*z);
    }
    void print()
    {
        qDebug()<<"("<<x<<","<<y<<","<<z<<")";
    }


};
struct Line{
    Point p1,p2; //有向线段，p1指向p2
    Line(){}
    Line(Point p1,Point p2){this->p1 = p1;this->p2 = p2; }
    void init(Point p1,Point p2){this->p1 = p1;this->p2 = p2; }
    bool operator ==(const Line &l)
    {
        return p1==l.p1&& p2 == l.p2;
    }
};


class Geometry
{
public:
    Geometry() {}
    static vec3 Cross(vec3 a, vec3 b)//法线
    {
        double x = a.y*b.z - a.z*b.y;
        double y = a.z*b.x - a.x*b.z;
        double z = a.x*b.y - a.y*b.x;
        return  vec3(x,y,z);
    }
    //计算距离
    static double PointToPoint(Point pos1, Point pos2)
    {
        double sum = 0;
        sum+= pow(pos2.x-pos1.x,2)+pow(pos2.y-pos1.y,2)+pow(pos2.z-pos1.z,2);
        return sqrt(sum);
    }
    static double PointToPoint(vec3 pos1, vec3 pos2)
    {
        double sum = 0;
        sum+= pow(pos2.x-pos1.x,2)+pow(pos2.y-pos1.y,2)+pow(pos2.z-pos1.z,2);
        return sqrt(sum);
    }
    //向量的模
    static double Norm(vec3 v){
        double sum = 0;
        sum+= pow(v.x,2)+pow(v.y,2)+pow(v.z,2);
        return sqrt(sum);
    }
    //点到线的距离
    static double PointToLine(Point p,Line l)
    {
        vec3 ab(l.p1,l.p2);
        vec3 ac(l.p1,p);
        return Norm(Cross(ab,ac))/PointToPoint(l.p1,l.p2); //d = (AB x AC)/|AB| ,|AB X AC|/2是三角形ABC的面积，这个三角形的底是|AB|，高就是C到AB的距离
    }

    static double PointToLine2D(Point p,Line l)
    {
        Point p1 = Point(l.p1.x, l.p1.y, 0);
        Point p2 = Point(l.p2.x, l.p2.y, 0);
        Point pp = Point(p.x, p.y, 0);
        vec3 ab(p1,p2);
        vec3 ac(p1,pp);
        return Norm(Cross(ab,ac))/PointToPoint(l.p1,l.p2); //d = (AB x AC)/|AB| ,|AB X AC|/2是三角形ABC的面积，这个三角形的底是|AB|，高就是C到AB的距离
    }

    //    static double angle(double x1,double y1, double x2,double y2){
    //        double n = (x1*x2+y1*y2);
    //        double m = sqrt(x1*x1+y1*y1)*sqrt(x2*x2+y2*y2);
    //        return acos(n/m)*180/M_PI;
    //    }
    static double angle3D(Point a,Point b,Point c){
        vec3 v1(a-b);
        vec3 v2(c-b);
        double x1 = v1.x;
        double y1 = v1.y;
        double z1 = v1.z;
        double x2 = v2.x;
        double y2 = v2.y;
        double z2 = v2.z;
        double n = (x1*x2 + y1*y2 + z1*z2);
        double m = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2 +z2*z2);
        return acos(n/m)*180/M_PI;
    }
    //空间两向量夹角
    static double angle3D(QVector3D n1,QVector3D n2)
    {
        double x1 = n1.x();
        double y1 = n1.y();
        double z1 = n1.z();
        double x2 = n2.x();
        double y2 = n2.y();
        double z2 = n2.z();
        double n = (x1*x2+y1*y2 + z1*z2);
        double m = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2 +z2*z2);
        return acos(n/m)*180/M_PI;
    }
    //空间两向量夹角
    static double angle3D(vec3 n1,vec3 n2)
    {
        double x1 = n1.x;
        double y1 = n1.y;
        double z1 = n1.z;
        double x2 = n2.x;
        double y2 = n2.y;
        double z2 = n2.z;
        double n = (x1*x2+y1*y2 + z1*z2);
        double m = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2 +z2*z2);
        return acos(n/m)*180/M_PI;
    }
    static double Cross2D(vec3 a, vec3 b)
    {
        return a.x*b.y-b.x*a.y;
    }
    static bool IsRightPoint(vec3 pt, Line line)
    {
        vec3 v1(line.p2.x-line.p1.x,line.p2.y-line.p1.y);//p1p2
        vec3 v2(line.p1.x-pt.x,line.p1.y-pt.y);// pp1
        double tmp =Cross2D(v1,v2);
        if(tmp>0) //大于0在右边
        {
            return true;
        }else
        {
            return false;
        }
    }
    static bool IsOnLine(vec3 pt, Line line)
    {
        vec3 v1(line.p2.x-line.p1.x,line.p2.y-line.p1.y);//p1p2
        vec3 v2(line.p1.x-pt.x,line.p1.y-pt.y);// pp1
        double tmp =Cross2D(v1,v2);
        double minx,miny;
        double maxx,maxy;
        if(line.p1.x < line.p2.x){
            minx = line.p1.x;
            maxx = line.p2.x;
        }
        else{
            minx = line.p2.x;
            maxx = line.p1.x;
        }

        if(line.p1.y < line.p2.y){
            miny= line.p1.y;
            maxy = line.p2.y;
        }
        else{
            miny = line.p2.y;
            maxy = line.p1.y;
        }
        if(fabs(tmp)<=1e-5 && pt.x>=minx&&pt.x <= maxx&& pt.y>=miny && pt.y <= maxy) //大于0在右边
        {
            return true;
        }else
        {
            return false;
        }
    }

    //    static vec3 calPointNormal(Point centerPoint, QList<Point> pts)
    //    {
    //        vec3 res;
    //        double tmp = (PointToPoint(pts[0], pts[1]))/(PointToPoint(centerPoint, pts[0]) + PointToPoint(centerPoint, pts[1]) + PointToPoint(pts[0], pts[1]));
    //        double f1, f2, f3;
    //        f1 = PointToPoint(pts[0], pts[1]);
    //        f2 = PointToPoint(centerPoint, pts[0]);
    //        f3 = PointToPoint(centerPoint, pts[1]);
    //        std::cout << f1 << "/(" << f2 << " + " << f3 << " + " << f1 << ")" << std::endl;
    //        std::cout << f1/(f1+f2+f3) << std::endl;
    //        res = Cross(vec3(pts[0], pts[1]), vec3(pts[0], pts[3])) * (PointToPoint(pts[0], pts[1])/(PointToPoint(centerPoint, pts[0]) + PointToPoint(centerPoint, pts[1]) + PointToPoint(pts[0], pts[1]))) +
    //                Cross(vec3(pts[1], centerPoint), vec3(pts[1], pts[3])) * (PointToPoint(pts[1], pts[2])/(PointToPoint(centerPoint, pts[1) + PointToPoint(centerPoint, pts[2]) + PointToPoint(pts[1], pts[2]))) +
    //                Cross(vec3(pts[0], pts[1]), vec3(pts[0], pts[3])) * (PointToPoint(pts[2], pts[3])/(PointToPoint(centerPoint, pts[0]) + PointToPoint(centerPoint, pts[1]) + PointToPoint(pts[0], pts[1])));
    //        qDebug() << tmp;
    //        return res;
    //    }
    static Point calCenterPoint(Point a, Point b, Point c)
    {
        return Point((a.x+b.x+c.x)/3, (a.y+b.y+c.y)/3, (a.z+b.z+c.z)/3);
    }
    static vec3 calNormalDiffAngle(Point target, QList<Point> pts, double *a1, double *a2)
    {
        QList<vec3> vecs;
        QList<vec3>::Iterator itor;
        //        vecs.push_back(Cross(vec3(pts[0], pts[1]), vec3(pts[0], pts[3])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[0], pts[1], pts[3]), target)));
        //        vecs.push_back(Cross(vec3(pts[1], target), vec3(pts[1], pts[3])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[3], pts[1], target), target)));
        //        vecs.push_back(Cross(vec3(pts[2], target), vec3(pts[2], pts[1])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[2], pts[1], target), target)));
        //        vecs.push_back(Cross(vec3(pts[4], target), vec3(pts[4], pts[2])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[2], pts[4], target), target)));
        //        vecs.push_back(Cross(vec3(pts[3], target), vec3(pts[3], pts[5])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[3], pts[5], target), target)));
        //        vecs.push_back(Cross(vec3(pts[5], target), vec3(pts[5], pts[6])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[5], pts[6], target), target)));
        //        vecs.push_back(Cross(vec3(pts[6], target), vec3(pts[6], pts[4])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[6], pts[4], target), target)));
        //        vecs.push_back(Cross(vec3(pts[7], pts[6]), vec3(pts[7], pts[4])) * (1.0/PointToPoint(Geometry::calCenterPoint(pts[7], pts[6], pts[4]), target)));
        vecs.push_back(Cross(vec3(pts[0], pts[1]), vec3(pts[0], pts[3])) * (PointToPoint(pts[0], target)/(PointToPoint(pts[0], target) + PointToPoint(pts[1], target) + PointToPoint(pts[0], pts[1]))));
        vecs.push_back(Cross(vec3(pts[1], target), vec3(pts[1], pts[3])) * (PointToPoint(pts[1], target)/(PointToPoint(pts[1], target) + PointToPoint(pts[3], target) + PointToPoint(pts[1], pts[3]))));
        vecs.push_back(Cross(vec3(pts[2], target), vec3(pts[2], pts[1])) * (PointToPoint(pts[2], target)/(PointToPoint(pts[2], target) + PointToPoint(pts[1], target) + PointToPoint(pts[2], pts[1]))));
        vecs.push_back(Cross(vec3(pts[4], target), vec3(pts[4], pts[2])) * (PointToPoint(pts[4], target)/(PointToPoint(pts[4], target) + PointToPoint(pts[2], target) + PointToPoint(pts[4], pts[2]))));
        vecs.push_back(Cross(vec3(pts[3], target), vec3(pts[3], pts[5])) * (PointToPoint(pts[3], target)/(PointToPoint(pts[3], target) + PointToPoint(pts[5], target) + PointToPoint(pts[3], pts[5]))));
        vecs.push_back(Cross(vec3(pts[5], target), vec3(pts[5], pts[6])) * (PointToPoint(pts[5], target)/(PointToPoint(pts[5], target) + PointToPoint(pts[6], target) + PointToPoint(pts[5], pts[6]))));
        vecs.push_back(Cross(vec3(pts[6], target), vec3(pts[6], pts[4])) * (PointToPoint(pts[6], target)/(PointToPoint(pts[6], target) + PointToPoint(pts[4], target) + PointToPoint(pts[6], pts[4]))));
        vecs.push_back(Cross(vec3(pts[7], pts[6]), vec3(pts[7], pts[4])) * (PointToPoint(pts[7], target)/(PointToPoint(pts[7], target) + PointToPoint(pts[4], target) + PointToPoint(pts[7], pts[4]))));
        itor = vecs.begin();
        vec3 res(0, 0, 0);
        while(itor != vecs.end())
        {
//            itor->print();
            res = res + *itor;
            itor++;
        }
        double theta, phi;
        theta = angle3D(vec3(1, 0, 0), vec3(res.x, res.y, 0));
        phi = angle3D(vec3(0, 0, 1), res);
//        res.print();
//        std::cout << theta << "------" << phi << std::endl;
        *a1 = theta;
        *a2 = phi;
        return res;
    }

    void static latLng2WebMercator(double lng, double lat, double* x, double* y)
    {
        double earthRad = 6378137.0;
        *x = (lng * M_PI * earthRad)/180;
        double a = lat * M_PI / 180;
        *y = earthRad / 2 * log((1.0 + sin(a)) / (1.0 - sin(a)));
    }

};
struct Triangle{
    Point p1,p2,p3;//三个点
    Point p[3];
    Line l1,l2,l3;//三条线
    Line l[3];//三条线
    Triangle(){}
    Triangle(Point a,Point b,Point c){ init(a,b,c);}
    bool isInTriangle(Point a)
    {
        bool r1 = Geometry::IsRightPoint(a,l1);
        bool r2 = Geometry::IsRightPoint(a,l2);
        bool r3 = Geometry::IsRightPoint(a,l3);
        if(a==p1 || a==p2 || a==p3)
            return false;
        if(r1 ==r2&& r2 == r3)
            return true;
        return false;
    }
    int isOnTriangle(Point a)
    {
        bool r1 = Geometry::IsOnLine(a,l1);
        bool r2 = Geometry::IsOnLine(a,l2);
        bool r3 = Geometry::IsOnLine(a,l3);
        if(a==p1 || a==p2 || a==p3)
            return false;
        if(r1)
            return 1;
        if(r2)
            return 2;
        if(r3)
            return 3;
        return 0;
    }

    void init(Point a,Point b,Point c){
        p1 =a ;
        p2 =b;
        p3 = c;
        p[0] = a;
        p[1] = b;
        p[2] = c;
        l1 = Line(a,b);
        l2 = Line(b,c);
        l3 = Line(c,a);
        l[0] = l1;
        l[1] = l2;
        l[2] = l3;
    }
    int  containsLine(Line l)
    {
        if((l.p1 == p1 && l.p2 == p2) || (l.p1 == p2 && l.p2 == p1) )
            return 1;
        if((l.p1 == p2 && l.p2 == p3) || (l.p1 == p3 && l.p2 == p2) )
            return 2;
        if((l.p1 == p1 && l.p2 == p3) || (l.p1 == p3 && l.p2 == p1) )
            return 3;
        return 0;
    }
    bool operator ==(const Triangle& tri)
    {
        if(p1 == tri.p1 && p2 == tri.p2&& p3 == tri.p3)
            return true;
        return false;
    }
    bool isTriangle()
    {
        if(p1==p2 || p1==p3 || p2==p3)
            return false;
        return true;
    }
};
struct Circle{
    double radius;
    vec3 center;
    Circle(){}
    Circle(vec3 cent,double r){ center = cent; radius = r;}
    static Circle genTriCircle(Triangle tri){
        Point p1 = tri.p1;
        Point p2 = tri.p2;
        Point p3 = tri.p3;
        double  x1,x2,x3,y1,y2,y3;
        x1  =  p1.x;
        x2  =  p2.x;
        x3  =  p3.x;
        y1  =  p1.y;
        y2  =  p2.y;
        y3  =  p3.y;
        //求外接圆半径
        double a=sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) );
        double b=sqrt( (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3) );
        double c=sqrt( (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3) );
        double p=(a+b+c)/2;
        double S=sqrt( p*(p-a)*(p-b)*(p-c) );

        //求外接圆圆心
        double t1=x1*x1+y1*y1;
        double t2=x2*x2+y2*y2;
        double t3=x3*x3+y3*y3;
        double temp=x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2;
        double x=(t2*y3+t1*y2+t3*y1-t2*y1-t3*y2-t1*y3)/temp/2;
        double y=(t3*x2+t2*x1+t1*x3-t1*x2-t2*x3-t3*x1)/temp/2;
        double r=a*b*c/(4*S);
        vec3 cent(x,y,0);
        return Circle(cent,r);

    }
    bool isInCircle(vec3 v)
    {
        double dist = Geometry::PointToPoint(center,v);
        return dist<radius;
    }
};



class Sort
{
public:
    Sort();

    static void quickSort(QVector<double>&a,int left,int right,bool(*cmp)(double ,double))
    {
        if(left<right){
            int i=left;
            int j=right;
            double  temp=a[left];
            if(left>=right)
                return;
            while(i<j)
            {
                while(i<j&& cmp(temp,a[j]))
                    j--;
                if(j>i)
                    a[i++]=a[j];//a[i]已经赋值给temp,所以直接将a[j]赋值给a[i],赋值完之后a[j],有空位
                while(i<j&&cmp(a[i],temp))
                    i++;
                if(i<j)
                    a[j--]=a[i];
            }
            a[i]=temp;//把基准插入,此时i与j已经相等R[low..pivotpos-1].keys≤R[pivotpos].key≤R[pivotpos+1..high].keys
            quickSort(a,left,i-1,cmp);/*递归左边*/
            quickSort(a,i+1,right,cmp);/*递归右边*/
        }
    }
    static void quickSort(QVector<Point>&a,int left,int right,bool(*cmp)(Point ,Point))
    {
        if(left<right){
            int i=left;
            int j=right;
            Point  temp=a[left];
            if(left>=right)
                return;
            while(i<j)
            {
                while(i<j&& cmp(temp,a[j]))
                    j--;
                if(j>i)
                    a[i++]=a[j];//a[i]已经赋值给temp,所以直接将a[j]赋值给a[i],赋值完之后a[j],有空位
                while(i<j&&cmp(a[i],temp))
                    i++;
                if(i<j)
                    a[j--]=a[i];
            }
            a[i]=temp;//把基准插入,此时i与j已经相等R[low..pivotpos-1].keys≤R[pivotpos].key≤R[pivotpos+1..high].keys
            quickSort(a,left,i-1,cmp);/*递归左边*/
            quickSort(a,i+1,right,cmp);/*递归右边*/
        }
    }
};
class Image{
public:
    static int getIndex(int i, int width)
    {
        if(i<0)
            i = 0;
        if(i>=width)
            i = width-1;
        return i;
    }
    static QImage TransToEdge(const QImage &source)
    {
        int w = source.width();
        int h = source.height();
        int Gmax =140;
        QImage Edge(w,h,QImage::Format_RGB32);

        for( int i = 0; i< h; i++){
            //卷积操作
            for(int j = 0; j < w; j++){
                double Gx =  (-1)* QColor(source.pixel(getIndex(j-1,w),getIndex(i-1,h))).red()
                        +(-2)*QColor(source.pixel(getIndex(j,w),getIndex(i-1,h))).red()
                        +(-1)*QColor(source.pixel(getIndex(j+1,w),getIndex(i-1,h))).red()
                        +QColor(source.pixel(getIndex(j-1,w),getIndex(i+1,h))).red()
                        +2*QColor(source.pixel(getIndex(j,w),getIndex(i+1,h))).red()
                        +QColor(source.pixel(getIndex(j+1,w),getIndex(i+1,h))).red();

                double Gy =  QColor(source.pixel(getIndex(j-1,w),getIndex(i-1,h))).red()
                        +(2)*QColor(source.pixel(getIndex(j-1,w),getIndex(i,h))).red()
                        +(1)*QColor(source.pixel(getIndex(j-1,w),getIndex(i+1,h))).red()
                        +(-1)*QColor(source.pixel(getIndex(j+1,w),getIndex(i-1,h))).red()
                        +(-2)*QColor(source.pixel(getIndex(j+1,w),getIndex(i,h))).red()
                        +(-1)*QColor(source.pixel(getIndex(j+1,w),getIndex(i+1,h))).red();

                double G = sqrt(Gx*Gx+Gy*Gy);

                QRgb pixel;
                if(G>Gmax)
                    pixel = qRgb(255,255,255);
                else
                    pixel = qRgb(0,0,0);
                Edge.setPixel(j,i,pixel);
            }
        }
        return Edge;
    }
};


#endif // UTIL_H
