#ifndef DEMREADER_H
#define DEMREADER_H
#include "opencv2/opencv.hpp"
#include "util.h"
#include <gdal_priv.h>
#define SIDE_THRESHOLD 5
#define THETA_THRESHOLD 20
#define PHI_THRESHOLD 15
#define MAX_COUNT_NUM 4
#define TILE_SIZE 256

typedef struct _TRIANGLE_DESC_
{
    cv::Point3d pt1, pt2, pt3;
    _TRIANGLE_DESC_(const cv::Point3d _pt1, const cv::Point3d _pt2, const cv::Point3d _pt3):
        pt1(_pt1), pt2(_pt2), pt3(_pt3){}
}TRIANGLE_DESC;

class DEMReader
{
public:
    DEMReader();
    ~DEMReader();
    std::vector<TRIANGLE_DESC> getCVTriangles(int startX, int startY, int sizeX, int sizeY, bool needFeature);
    void getGoogleMapTile(double lng, double lat, cv::Point *pixelCoord, cv::Point *tileCoord);
    void getGoogleMapTile(double lng, double lat, cv::Point *tileCoord);
    void calTextureRange();
    int getWidth();
    int getHeight();
    double getOriginX();
    double getOriginY();
    double getDiffX();
    double getDiffY();

private:
    cv::Point2d project(double lat, double lng);
    QList<Point> getPoints(int startX, int startY, int sizeX, int sizeY);
    QList<Point> featurePointSelection(QList<Point> pts, int row, int col);
    std::vector<TRIANGLE_DESC> delaunayAlgorithm(const cv::Rect boundRc,const std::vector<cv::Point3d>& points, QVector<double> elev, int col, int xOff, int yOff);
    double originX;//dem data start X
    double originY;//dem data start Y
    double diffX;//x intervals between pixels
    double diffY;//y intervals between pixels
    double startX;//my dem reading start from x
    double startY;//my dem reading start from y
    double sizeX;//my dem reading scale x
    double sizeY;//my dem reading scale y
    int height;
    int width;
    int texStartXNum;
    int texStartYNum;
    int texEndXNum;
    int texEndYNum;
    GDALRasterBand *poBand;
    GDALDataset *poDataset;

};

#endif // DEMREADER_H
