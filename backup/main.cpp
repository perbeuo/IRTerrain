#include <iostream>
#include <stdint.h>
#include <math.h>
#include <iomanip>
#include <osgViewer/Viewer>
#include <osg/Geode>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osg/Texture2D>
#include <osg/TexGen>
#include <osg/ShapeDrawable>
#include <gdal_priv.h>
#include <osgTerrain/TerrainTile>
#include <osgTerrain/GeometryTechnique>
#include <osgTerrain/Layer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Shape>
#include <osg/PolygonMode>
#include <osgUtil/DelaunayTriangulator>
#include <osgUtil/Optimizer>
#include <osgUtil/SmoothingVisitor>
#include <osg/Geometry>
#include "opencv2/opencv.hpp"
#include "util.h"
#define SIDE_THRESHOLD 5
#define THETA_THRESHOLD 20
#define PHI_THRESHOLD 15
#define MAX_COUNT_NUM 4
#define TILE_SIZE 256
using namespace std;
typedef struct _TRIANGLE_DESC_
{
    cv::Point3d pt1, pt2, pt3;
    //    vec3 norm1, norm2, norm3;
    _TRIANGLE_DESC_(const cv::Point3d _pt1, const cv::Point3d _pt2, const cv::Point3d _pt3/*, const vec3 _norm1, const vec3 _norm2, const vec3 _norm3*/):
        pt1(_pt1), pt2(_pt2), pt3(_pt3)/*, norm1(_norm1), norm2(_norm2), norm3(_norm3)*/{}
}TRIANGLE_DESC;

vector<TRIANGLE_DESC> delaunayAlgorithm(const cv::Rect boundRc, const vector<cv::Point3d>& points, QVector<double> elev, /*QList<vec3> norms,*/ int col, int xOff, int yOff);
QList<Point> featurePointSelection(QList<Point> pts, int row, int col/*, QList<vec3> *resNorms*/);

cv::Point2d project(double lat, double lng) {
    double siny = sin(lat * M_PI / 180);

    // Truncating to 0.9999 effectively limits latitude to 89.189. This is
    // about a third of a tile past the edge of the world tile.
    siny = min(max(siny, -0.9999), 0.9999);
    double x = TILE_SIZE * (0.5 + lng / 360);
    double y = TILE_SIZE * (0.5 - log((1 + siny) / (1 - siny)) / (4 * M_PI));

    return cv::Point2d(x, y);
}

QList<Point> getPoints(int startX, int startY, int sizeX, int sizeY)
{
    //    const char *pszFilename="/home/lzt/apps/tin-terrain/3rdparty/craterlake/dems_10m.dem";
    const char *pszFilename="/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif";
    //"/home/lzt/material/DEM/n39_e126_1arc_v3.tif";
    double originX, originY, diffX, diffY;
    GDALDataset *poDataset;
    QList<Point> pts;
    GDALAllRegister();
    poDataset = (GDALDataset*)GDALOpen(pszFilename, GA_ReadOnly);
    if (poDataset == NULL)
    {
        cout << "open failed!" << endl;
    }

    double adfGeoTransform[6];
    printf("Driver: %s/%s\n",
           poDataset->GetDriver()->GetDescription(),
           poDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
    printf("Size is: %dx%dx%d\n",
           poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
           poDataset->GetRefCount());
    if (poDataset->GetProjectionRef() != NULL)
        printf("Projection is `%s`\n", poDataset->GetProjectionRef());
    if(poDataset->GetGeoTransform(adfGeoTransform) == CE_None)
    {
        printf("Origin = (%.15f, %.15f)\n",
               adfGeoTransform[0], adfGeoTransform[3]);
        printf("Pixel Size = (%.15f, %.15f)\n",
               adfGeoTransform[1], adfGeoTransform[5]);
        //        printf("2 4 = (%.15f, %.15f)\n",
        //               adfGeoTransform[2], adfGeoTransform[4]);
    }
    originX = adfGeoTransform[0];
    originY = adfGeoTransform[3];
    diffX = adfGeoTransform[1];
    diffY = adfGeoTransform[5];

    GDALRasterBand *poBand;
    int nBlockXSize, nBlockYSize;
    int bGotMin, bGotMax;
    double adfMinMax[2];
    poBand = poDataset->GetRasterBand(1);
    poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
    printf("Block=%dx%d  Type=%s  ColorInterp=%s\n",
           nBlockXSize, nBlockYSize,
           GDALGetDataTypeName(poBand->GetRasterDataType()),
           GDALGetColorInterpretationName(poBand->GetColorInterpretation()));
    adfMinMax[0] = poBand->GetMinimum(&bGotMin);
    adfMinMax[1] = poBand->GetMaximum(&bGotMax);
    if(!(bGotMin && bGotMax))
        GDALComputeRasterMinMax((GDALRasterBandH)poBand,
                                TRUE, adfMinMax);
    printf("Min=%.3f, Max=%.3f\n", adfMinMax[0],
            adfMinMax[1]);
    if(poBand->GetOverviewCount() > 0)
        printf("Band has %d overviews.\n",
               poBand->GetOverviewCount());
    if(poBand->GetColorTable() != NULL)
        printf("Band has a color table with %d entries.\n",
               poBand->GetColorTable()->GetColorEntryCount());

    float *pafScanline;
    //    int nXSize = poBand->GetXSize();
    //    int nYSize = poBand->GetYSize();
    int nXSize = sizeX;
    int nYSize = sizeY;
    pafScanline = (float *)CPLMalloc(sizeof(float)*nXSize*nYSize);
    poBand->RasterIO(GF_Read, startX, startY, nXSize, nYSize,
                     pafScanline, nXSize, nYSize, GDT_Float32,
                     0, 0);
    for (int i = 0; i < nYSize; i++)
    {
        for (int j = 0; j < nXSize; j++)
        {
            double x, y;
            double elevation = *pafScanline;
            //            printf("%.3f\t", elevation);
            //            pts.push_back(Point(startX+i*diffX, startY+j*diffY, elevation));
            Geometry::latLng2WebMercator(adfGeoTransform[0]+(j+startX)*adfGeoTransform[1]+(i+startY)*adfGeoTransform[2],
                    adfGeoTransform[3] + (j+startX) * adfGeoTransform[4] + (i+startY) * adfGeoTransform[5], &x, &y);
            //            x = adfGeoTransform[0]+(i+startX)*adfGeoTransform[1]+(j+startY)*adfGeoTransform[2];
            //            y = adfGeoTransform[3] + (i+startX) * adfGeoTransform[4] + (j+startY) * adfGeoTransform[5];
            pts.push_back(Point(x, y, elevation));
            pafScanline++;
        }
        //        printf("\n");
    }
    //    double x1, y1, x2, y2;
    //    Geometry::latLng2WebMercator( startX,  startY, &x1, &y1);
    //    Geometry::latLng2WebMercator( startX +  diffX,  startY, &x2, &y2);
    //    cout << setprecision(15) << x2-x1 << "," << y2-y1 << endl;
    //    cout << pts.size() << endl;
    //    pts[3600].print();
    //    pts[3601].print();
    //    pts[3602].print();

    QList<Point> adjusted;
    for (int i = nYSize-1; i >= 0; i--)
    {
        for (int j = 0; j < nXSize; j++)
        {
            adjusted.push_back(pts[j*nYSize+i]);
        }
    }
    //    for (int i = 0; i < 100; i++)
    //    {
    //        pts[i].print();
    //    }
    GDALClose((GDALDatasetH)poDataset);

    return pts;
}

vector<TRIANGLE_DESC> getCVTriangles(int startX, int startY, int sizeX, int sizeY/*, QList<vec3> *resNorms*/, bool needFeature)
{
    const char *pszFilename="/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif";
    double originX, originY, diffX, diffY;
    GDALDataset *poDataset;
    QList<Point> pts;
    QVector<double> elev;
    GDALAllRegister();
    poDataset = (GDALDataset*)GDALOpen(pszFilename, GA_ReadOnly);
    if (poDataset == NULL)
    {
        cout << "open failed!" << endl;
    }

    double adfGeoTransform[6];
    printf("Driver: %s/%s\n",
           poDataset->GetDriver()->GetDescription(),
           poDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
    printf("Size is: %dx%dx%d\n",
           poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
           poDataset->GetRefCount());
    if (poDataset->GetProjectionRef() != NULL)
        printf("Projection is `%s`\n", poDataset->GetProjectionRef());
    if(poDataset->GetGeoTransform(adfGeoTransform) == CE_None)
    {
        printf("Origin = (%.15f, %.15f)\n",
               adfGeoTransform[0], adfGeoTransform[3]);
        printf("Pixel Size = (%.15f, %.15f)\n",
               adfGeoTransform[1], adfGeoTransform[5]);
        //        printf("2 4 = (%.15f, %.15f)\n",
        //               adfGeoTransform[2], adfGeoTransform[4]);
    }
    originX = adfGeoTransform[0];
    originY = adfGeoTransform[3];
    diffX = adfGeoTransform[1];
    diffY = adfGeoTransform[5];

    GDALRasterBand *poBand;
    int nBlockXSize, nBlockYSize;
    int bGotMin, bGotMax;
    double adfMinMax[2];
    poBand = poDataset->GetRasterBand(1);
    poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
    printf("Block=%dx%d  Type=%s  ColorInterp=%s\n",
           nBlockXSize, nBlockYSize,
           GDALGetDataTypeName(poBand->GetRasterDataType()),
           GDALGetColorInterpretationName(poBand->GetColorInterpretation()));
    adfMinMax[0] = poBand->GetMinimum(&bGotMin);
    adfMinMax[1] = poBand->GetMaximum(&bGotMax);
    if(!(bGotMin && bGotMax))
        GDALComputeRasterMinMax((GDALRasterBandH)poBand,
                                TRUE, adfMinMax);
    printf("Min=%.3f, Max=%.3f\n", adfMinMax[0],
            adfMinMax[1]);
    if(poBand->GetOverviewCount() > 0)
        printf("Band has %d overviews.\n",
               poBand->GetOverviewCount());
    if(poBand->GetColorTable() != NULL)
        printf("Band has a color table with %d entries.\n",
               poBand->GetColorTable()->GetColorEntryCount());

    float *pafScanline;
    //    int nXSize = poBand->GetXSize();
    //    int nYSize = poBand->GetYSize();
    int nXSize = sizeX;
    int nYSize = sizeY;
    pafScanline = (float *)CPLMalloc(sizeof(float)*nXSize*nYSize);
    poBand->RasterIO(GF_Read, startX, startY, nXSize, nYSize,
                     pafScanline, nXSize, nYSize, GDT_Float32,
                     0, 0);
    for (int i = 0; i < nYSize; i++)
    {
        for (int j = 0; j < nXSize; j++)
        {
            double x, y;
            double elevation = *pafScanline;
            //            printf("%.3f\t", elevation);
            //            pts.push_back(Point(startX+i*diffX, startY+j*diffY, elevation));
            //            Geometry::latLng2WebMercator(adfGeoTransform[0]+(j+startX)*adfGeoTransform[1]+(i+startY)*adfGeoTransform[2],
            //                    adfGeoTransform[3] + (j+startX) * adfGeoTransform[4] + (i+startY) * adfGeoTransform[5], &x, &y);
            //            x = adfGeoTransform[0]+(j+startX)*adfGeoTransform[1]+(i+startY)*adfGeoTransform[2];
            //            y = adfGeoTransform[3] + (j+startX) * adfGeoTransform[4] + (i+startY) * adfGeoTransform[5];
            x = (j+startX)*30;
//            y = (nYSize+startY-1-i)*30;
            y = (i+startY)*30;
            pts.push_back(Point(x, y, elevation));
            elev.push_back(elevation);
            pafScanline++;
        }
        //        printf("\n");
    }
    //    double x1, y1, x2, y2;
    //    Geometry::latLng2WebMercator( startX,  startY, &x1, &y1);
    //    Geometry::latLng2WebMercator( startX +  diffX,  startY, &x2, &y2);
    //    cout << setprecision(15) << x2-x1 << "," << y2-y1 << endl;
    //    cout << pts.size() << endl;
    //    pts[3600].print();
    //    pts[3601].print();
    //    pts[3602].print();

    //    QList<Point> adjusted;
    //    for (int i = nYSize-1; i >= 0; i--)
    //    {
    //        for (int j = 0; j < nXSize; j++)
    //        {
    //            adjusted.push_back(pts[j*nYSize+i]);
    //        }
    //    }
    //    for (int i = 0; i < 100; i++)
    //    {
    //        pts[i].print();
    //    }
    //Google map tile test
    double lng = adfGeoTransform[0]+(startX)*adfGeoTransform[1]+(startY)*adfGeoTransform[2];
    double lat = adfGeoTransform[3]+(startX)*adfGeoTransform[4]+(startY)*adfGeoTransform[5];
    cv::Point2d worldCoord = project(lat, lng);
    cout << "world:" << worldCoord.x << ", " << worldCoord.y << endl;
    int scale = 1 << 10;
    cv::Point pixelCoord = cv::Point(floor(worldCoord.x * scale), floor(worldCoord.y * scale));
    cout << "pixel:" << pixelCoord.x << ", " << pixelCoord.y << endl;
    cv::Point tileCoord = cv::Point(floor(worldCoord.x * scale/TILE_SIZE), floor(worldCoord.y * scale/TILE_SIZE));
    cout << "tile:" << tileCoord.x << ", " << tileCoord.y << endl;
    cout << "pixel in image:" << pixelCoord.x-tileCoord.x*TILE_SIZE << ", " << pixelCoord.y-tileCoord.y*TILE_SIZE << endl;

    lng = adfGeoTransform[0]+(startX+sizeX)*adfGeoTransform[1]+(startY+sizeY)*adfGeoTransform[2];
    lat = adfGeoTransform[3]+(startX+sizeX)*adfGeoTransform[4]+(startY+sizeY)*adfGeoTransform[5];
    worldCoord = project(lat, lng);
    cout << "world:" << worldCoord.x << ", " << worldCoord.y << endl;
    pixelCoord = cv::Point(floor(worldCoord.x * scale), floor(worldCoord.y * scale));
    cout << "pixel:" << pixelCoord.x << ", " << pixelCoord.y << endl;
    tileCoord = cv::Point(floor(worldCoord.x * scale/TILE_SIZE), floor(worldCoord.y * scale/TILE_SIZE));
    cout << "tile:" << tileCoord.x << ", " << tileCoord.y << endl;
    cout << "pixel in image:" << pixelCoord.x-tileCoord.x*TILE_SIZE << ", " << pixelCoord.y-tileCoord.y*TILE_SIZE << endl;

    GDALClose((GDALDatasetH)poDataset);
    QList<Point> feature = featurePointSelection(pts, sizeX, sizeY/*, resNorms*/);

    //opencv delaunay
    const int width = sizeX*32;
    const int height = sizeY*32;
    double offX, offY;
    //    Geometry::latLng2WebMercator(adfGeoTransform[0]+(startX)*adfGeoTransform[1]+(startY)*adfGeoTransform[2],
    //            adfGeoTransform[3] + (startX) * adfGeoTransform[4] + (startY) * adfGeoTransform[5], &offX, &offY);
    offX = startX*30;
    offY = startY*30;
    vector<cv::Point3d> testPoints;
    if (needFeature)
    {
        for (int i = 0; i < feature.size(); i++)
        {
            testPoints.push_back(cv::Point3d(feature[i].x, feature[i].y, feature[i].z));
        }
    }else
    {
        for (int i = 0; i < pts.size(); i++)
        {
            testPoints.push_back(cv::Point3d(pts[i].x, pts[i].y, pts[i].z));
        }
    }
    //    const cv::Rect pageRc(offX-1, offY-height+1, width, height);
    const cv::Rect pageRc(offX-1, offY-1, width, height);
    const auto triangles = delaunayAlgorithm(pageRc,testPoints, elev,/* *resNorms,*/ nXSize, offX, offY);

    return triangles;
}

osg::Node* createHeightField(std::string heightFile, std::string texFile) {

    osg::Image* heightMap = osgDB::readImageFile(heightFile);

    osg::HeightField* heightField = new osg::HeightField();
    heightField->allocate(heightMap->s(), heightMap->t());
    heightField->setOrigin(osg::Vec3(-heightMap->s() / 2, -heightMap->t() / 2, 0));
    heightField->setXInterval(1.0f);
    heightField->setYInterval(1.0f);
    heightField->setSkirtHeight(1.0f);

    for (int r = 0; r < heightField->getNumRows(); r++) {
        for (int c = 0; c < heightField->getNumColumns(); c++) {
            heightField->setHeight(c, r, ((*heightMap->data(c, r)) / 255.0f) * 80.0f);
        }
    }

    osg::Geode* geode = new osg::Geode();
    geode->addDrawable(new osg::ShapeDrawable(heightField));

    //    osg::Texture2D* tex = new osg::Texture2D(osgDB::readImageFile(texFile));
    //    tex->setFilter(osg::Texture2D::MIN_FILTER,osg::Texture2D::LINEAR_MIPMAP_LINEAR);
    //    tex->setFilter(osg::Texture2D::MAG_FILTER,osg::Texture2D::LINEAR);
    //    tex->setWrap(osg::Texture::WRAP_S, osg::Texture::REPEAT);
    //    tex->setWrap(osg::Texture::WRAP_T, osg::Texture::REPEAT);
    //    geode->getOrCreateStateSet()->setTextureAttributeAndModes(0, tex);

    return geode;
}

osg::ref_ptr<osgTerrain::TerrainTile> readFromDEM()
{
    GDALAllRegister();
    GDALDataset* poDataset=(GDALDataset*)GDALOpen("/home/lzt/material/DEM/crater_dems_10m.dem",GA_ReadOnly);
    if(poDataset){
        double gdalGeoTransform[6];
        poDataset->GetGeoTransform(gdalGeoTransform);
        osg::HeightField* hf=new osg::HeightField();
        hf->allocate(poDataset->GetRasterXSize(),poDataset->GetRasterYSize());
        hf->setOrigin(osg::Vec3(gdalGeoTransform[0],gdalGeoTransform[3],0));
        hf->setXInterval(gdalGeoTransform[2]);
        hf->setYInterval(gdalGeoTransform[5]);
        float * heightData=new float[poDataset->GetRasterXSize()*poDataset->GetRasterYSize()];
        poDataset->GetRasterBand(1)->RasterIO(GF_Read,0,0,poDataset->GetRasterXSize(),poDataset->GetRasterYSize(),heightData,poDataset->GetRasterXSize(),poDataset->GetRasterYSize(),GDT_Float32,0,0);
        float* heightPtr=heightData;
        float noDataValueFill=0.0f;
        float noDataValue=poDataset->GetRasterBand(1)->GetNoDataValue();
        for(int r=poDataset->GetRasterYSize()-1;r>=0;--r){
            for(int c=0;c<poDataset->GetRasterXSize();++c){
                float h=*heightPtr++;
                if(h!=noDataValue)
                    hf->setHeight(c,r,h);
                else
                    hf->setHeight(c,r,noDataValueFill);
            }
        }

        osg::ref_ptr<osgTerrain::TerrainTile> terrainTile=new osgTerrain::TerrainTile;
        osg::ref_ptr<osgTerrain::Locator> locator=new osgTerrain::Locator;
        double minX,minY,maxX,maxY;
        minX=std::min(gdalGeoTransform[0],gdalGeoTransform[0]+poDataset->GetRasterXSize()*gdalGeoTransform[1]);
        minY=std::min(gdalGeoTransform[3],gdalGeoTransform[3]+poDataset->GetRasterYSize()*gdalGeoTransform[5]);
        maxX=std::max(gdalGeoTransform[0],gdalGeoTransform[0]+poDataset->GetRasterXSize()*gdalGeoTransform[1]);
        maxY=std::max(gdalGeoTransform[3],gdalGeoTransform[3]+poDataset->GetRasterYSize()*gdalGeoTransform[5]);
        locator->setTransformAsExtents( minX,minY,maxX,maxY);

        osg::ref_ptr<osgTerrain::HeightFieldLayer> hfl=new osgTerrain::HeightFieldLayer;
        hfl->setHeightField(hf);
        hfl->setLocator(locator.get());
        terrainTile->setElevationLayer(hfl);
        osg::ref_ptr<osg::PolygonMode> pLolyMode= new osg::PolygonMode;
        osg::ref_ptr<osg::StateSet> stateset = new osg::StateSet;
        pLolyMode->setMode(osg::PolygonMode::FRONT_AND_BACK ,osg::PolygonMode::LINE);
        stateset->setAttribute(pLolyMode);
        terrainTile->setStateSet(stateset);
        return terrainTile;
    }
    return NULL;
}

osg::Drawable *createTriangulate(QList<Point> pts)
{
    //    osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array();
    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array();
    //    osg::ref_ptr<osg::Vec3Array> norms = new osg::Vec3Array();
    //    float vertex[][3] = { -5.0f,-5.0f, 0.4f,
    //                          1.0f, -5.6f, 0.0f,
    //                          5.0f, -4.0f, -0.5f,
    //                          -6.2f, 0.0f, 4.2f,
    //                          -1.0f,-0.5f, 4.8f,
    //                          4.3f, 1.0f, 3.0f,
    //                          -4.8f, 5.4f, 0.3f,
    //                          0.6f, 5.1f,-0.8f,
    //                          5.2f, 4.5f, 0.1f };
    //    unsigned int n = sizeof(vertex) / sizeof(float[3]);
    for (unsigned int i = 0; i < pts.size(); i++)
    {
        //        coords->push_back(osg::Vec3(pts[i].x, pts[i].y, 0));
        points->push_back(osg::Vec3(pts[i].x, pts[i].y, pts[i].z));
    }
    cout << pts.size() << "==" << points->size() << endl;
    cout << "converting..." << endl;
    osg::ref_ptr<osgUtil::DelaunayTriangulator> dt = new osgUtil::DelaunayTriangulator();
    dt->setInputPointArray(points);
    //    dt->setOutputNormalArray(norms);
    dt->triangulate();
    cout << "converting complete" << endl;
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
    geometry->setVertexArray(points.get());
    //    geometry->setNormalArray(norms.get());
    //    geometry->setNormalBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);

    geometry->addPrimitiveSet(dt->getTriangles());
    //    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points->size()));
    return geometry.release();

}

QList<Point> featurePointSelection(QList<Point> pts, int row, int col/*, QList<vec3> *resNorms*/)
{
    QList<Point> respts;
    vec3 norms[row][col];
    //    QList<vec3> resNorms;
    double thetas[row][col];
    double phis[row][col];

    //    resNorms->clear();
    //cal normal vectors
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if(i == 0 || i == row-1 || j == 0 || j == col-1)
            {
                thetas[i][j] = DBL_MAX;
                phis[i][j] = DBL_MAX;
            }
            else
            {
                double theta, phi;
                QList<Point> surroundPts;
                surroundPts.push_back(pts[(i-1)*col+j-1]);
                surroundPts.push_back(pts[(i-1)*col+j]);
                surroundPts.push_back(pts[(i-1)*col+j+1]);
                surroundPts.push_back(pts[i*col+j-1]);
                surroundPts.push_back(pts[i*col+j+1]);
                surroundPts.push_back(pts[(i+1)*col+j-1]);
                surroundPts.push_back(pts[(i+1)*col+j]);
                surroundPts.push_back(pts[(i+1)*col+j+1]);
                vec3 norm = Geometry::calNormalDiffAngle(pts[i*col+j], surroundPts, &theta, &phi);
                norms[i][j] = norm;
                thetas[i][j] = theta;
                phis[i][j] = phi;
            }
        }
    }
    //    for (int i = 0; i < row; i++)
    //    {
    //        for (int j = 0; j < col; j++)
    //        {
    //            std::cout << thetas[i][j] << "   ";
    //        }
    //        std::cout << std::endl;
    //    }

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            //for the first and last row
            if (i == 0 || i == row-1)
            {
                //add corner
                if (j == 0 || j == col-1)
                {
                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
                }
                else
                {
                    double diff1, diff2, diff3;
                    //same or different direction of slope means different gradient
                    diff1 = pts[i*col+j].z - pts[i*col+j-1].z;
                    diff2 = pts[i*col+j+1].z - pts[i*col+j].z;
                    diff3 = fabs(diff1-diff2);
                    //not the corner point
                    if (diff3 > SIDE_THRESHOLD)
                    {
                        respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
                    }

                }
                //                resNorms->push_back(vec3(0, 0, 0));
            }
            //for the first and last column(without cornor points)
            else if (j == 0 || j == col-1)
            {
                double diff1, diff2, diff3;
                //same or different direction of slope means different gradient
                diff1 = pts[i*col+j].z - pts[(i-1)*col+j].z;
                diff2 = pts[(i+1)*col+j].z - pts[i*col+j].z;
                diff3 = fabs(diff1-diff2);
                if (diff3 > SIDE_THRESHOLD)
                {
                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
                }
                //                resNorms->push_back(vec3(0, 0, 0));
            }
            //for other points
            else
            {
                QList<double> diffThetas, diffPhis, tmpThetas, tmpPhis;

                tmpThetas.push_back(thetas[i-1][j-1]);
                tmpThetas.push_back(thetas[i-1][j]);
                tmpThetas.push_back(thetas[i-1][j+1]);
                tmpThetas.push_back(thetas[i][j-1]);
                tmpThetas.push_back(thetas[i][j+1]);
                tmpThetas.push_back(thetas[i+1][j-1]);
                tmpThetas.push_back(thetas[i+1][j]);
                tmpThetas.push_back(thetas[i+1][j+1]);

                tmpPhis.push_back(phis[i-1][j-1]);
                tmpPhis.push_back(phis[i-1][j]);
                tmpPhis.push_back(phis[i-1][j+1]);
                tmpPhis.push_back(phis[i][j-1]);
                tmpPhis.push_back(phis[i][j+1]);
                tmpPhis.push_back(phis[i+1][j-1]);
                tmpPhis.push_back(phis[i+1][j]);
                tmpPhis.push_back(phis[i+1][j+1]);

                if (i == 4 && j == 134)
                {
                    cout << "here" << endl;
                }
                for (int k = 0; k < tmpThetas.size(); k++)
                {
                    double difftheta, diffphi;

                    if (fabs(tmpThetas[k] - DBL_MAX) > 1e-6)
                    {
                        difftheta = fabs(thetas[i][j] - tmpThetas[k]);
                        if (difftheta > 180)
                            difftheta = 360 - difftheta;
                        diffThetas.push_back(difftheta);
                    }
                    if (fabs(tmpPhis[k] - DBL_MAX) > 1e-6)
                    {
                        diffphi = fabs(phis[i][j] - tmpPhis[k]);
                        diffPhis.push_back(diffphi);
                    }
                }
                //                cout << "Point at (" << i << "," << j << ")" << endl;
                int count = 0;
                for(int k = 0; k < diffThetas.size(); k++)
                {
                    //                    cout << "theta:" << diffThetas[k] << endl;
                    //                    cout << "phi:" << diffPhis[k] << endl;
                    if (diffThetas[k] > THETA_THRESHOLD)
                        count++;
                    if (diffPhis[k] > PHI_THRESHOLD)
                        count++;
                }
                if (count > MAX_COUNT_NUM)
                {
                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
                }
                //                resNorms->push_back(norms[i][j]);
            }
        }
    }
    cout << "feature point num:" << respts.size() << endl;
    //    for (int i = 0; i < respts.size(); i++)
    //        cout<<"("<<respts[i].x<<","<<respts[i].y<<","<<respts[i].z<<")"<<endl;
    return respts;
}

vector<TRIANGLE_DESC> delaunayAlgorithm(const cv::Rect boundRc,const vector<cv::Point3d>& points, QVector<double> elev, /*QList<vec3> norms, */int col, int xOff, int yOff)
{

    if (points.empty())
    {
        return vector<TRIANGLE_DESC>();
    }
    vector<TRIANGLE_DESC> result;

    vector<cv::Vec6f> _temp_result;
    cv::Subdiv2D subdiv2d(boundRc);
    //    double dis=1e-6;
    for (const auto point : points)
    {
        subdiv2d.insert(cv::Point2d(point.x, point.y));
    }
    subdiv2d.getTriangleList(_temp_result);

    for (const auto _tmp_vec : _temp_result)
    {
        double z1, z2, z3;
        //        int count = 0;
        //        for(int i = 0; i < points.size(); i++)
        //        {
        //            if(fabs(points.at(i).x-_tmp_vec[0]) < dis &&
        //                    fabs(points.at(i).y-_tmp_vec[1]) < dis)
        //            {
        //                z1 = points.at(i).z;
        //                count++;
        //            }else if(fabs(points.at(i).x-_tmp_vec[2]) < dis &&
        //                     fabs(points.at(i).y-_tmp_vec[3]) < dis)
        //            {
        //                z2 = points.at(i).z;
        //                count++;
        //            }else if(fabs(points.at(i).x-_tmp_vec[4]) < dis &&
        //                     fabs(points.at(i).y-_tmp_vec[5]) < dis)
        //            {
        //                z3 = points.at(i).z;
        //                count++;
        //            }
        //            if (count == 3)
        //                break;
        //        }
        z1 = elev[(_tmp_vec[0]-xOff)/30 + ((_tmp_vec[1]-yOff)/30)*col];
        z2 = elev[(_tmp_vec[2]-xOff)/30 + ((_tmp_vec[3]-yOff)/30)*col];
        z3 = elev[(_tmp_vec[4]-xOff)/30 + ((_tmp_vec[5]-yOff)/30)*col];
        cv::Point3d pt1(_tmp_vec[0], _tmp_vec[1], z1);
        cv::Point3d pt2(_tmp_vec[2], _tmp_vec[3], z2);
        cv::Point3d pt3(_tmp_vec[4], _tmp_vec[5], z3);
        //        vec3 norm1 = norms[(_tmp_vec[0]-xOff)/30 + ((_tmp_vec[1]-yOff)/30)*col];
        //        vec3 norm2 = norms[(_tmp_vec[2]-xOff)/30 + ((_tmp_vec[3]-yOff)/30)*col];
        //        vec3 norm3 = norms[(_tmp_vec[2]-xOff)/30 + ((_tmp_vec[3]-yOff)/30)*col];
        result.push_back(TRIANGLE_DESC(pt1, pt2, pt3/*, norm1, norm2, norm3*/));
    }
    return result;
}

//创建二维纹理状态对象
osg::ref_ptr<osg::Texture2D> createTexture2D(osg::ref_ptr<osg::Image> image)
{
    //创建二维纹理对象
    osg::ref_ptr<osg::Texture2D> texture = new osg::Texture2D();
    texture->setDataVariance(osg::Object::DYNAMIC);
    //设置贴图
    texture->setImage(image.get());
    //设置滤波
    texture->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);
    texture->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
    //设置环绕模式
    texture->setWrap(osg::Texture::WRAP_S, osg::Texture::REPEAT);
    texture->setWrap(osg::Texture::WRAP_T, osg::Texture::REPEAT);

    //创建状态集对象
//    osg::ref_ptr<osg::StateSet> stateset = new osg::StateSet();
//    stateset->setTextureAttributeAndModes(0, texture.get(), osg::StateAttribute::ON);

    return texture.release();
}

osg::ref_ptr<osg::Geode> createTerrain()
{
    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
    osg::ref_ptr<osg::Image> image = osgDB::readImageFile("/home/lzt/material/img_map/gs_637_474_10.jpg");
    int xSize = 100, ySize = 100;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if (i==1&&j==1)
                break;
            vector<TRIANGLE_DESC> triangles;
            osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
            triangles = getCVTriangles(i*xSize, j*ySize, xSize+1, ySize+1, false);
            osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array();
            osg::ref_ptr<osg::Vec2Array> texCoord = new osg::Vec2Array();

            for (const auto triangle : triangles)
            {
                points->push_back(osg::Vec3(triangle.pt1.x, triangle.pt1.y, triangle.pt1.z));
                points->push_back(osg::Vec3(triangle.pt2.x, triangle.pt2.y, triangle.pt2.z));
                points->push_back(osg::Vec3(triangle.pt3.x, triangle.pt3.y, triangle.pt3.z));
//                double scale = triangle.pt3.x/(30*xSize);
                texCoord->push_back(osg::Vec2(triangle.pt1.x/(30*(xSize+1)), triangle.pt1.y/(30*(ySize+1))));
                texCoord->push_back(osg::Vec2(triangle.pt2.x/(30*(xSize+1)), triangle.pt2.y/(30*(ySize+1))));
                texCoord->push_back(osg::Vec2(triangle.pt3.x/(30*(xSize+1)), triangle.pt3.y/(30*(ySize+1))));
            }

            geometry->setVertexArray(points.get());
            geometry->setTexCoordArray(0,texCoord.get());
            cout << "pt size" << points->size() << endl;
            geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points->size()));
            geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0, createTexture2D(image), osg::StateAttribute::ON);

            osgUtil::SmoothingVisitor smv;
            smv.smooth(*geometry);

            geode->addDrawable(geometry.get());
        }
    }
    return geode.release();
}

int main()
{
    //    osgViewer::Viewer viewer;
    //    osg::Group *root = new osg::Group();
    //    //    osg::ref_ptr <osg::Node> b25 = createHeightField("/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif", "");
    //    osg::ref_ptr <osg::Node> b25 = createHeightField("/home/lzt/material/DEM/n41_e129_1arc_v3.tif", "");
    //    root->addChild(b25.get());
    //    viewer.setSceneData(root);
    //    viewer.realize();
    //    return viewer.run();

    //    QList<Point> pts = getPoints(0, 0, 100, 100);
    //    QList<Point> feature = featurePointSelection(pts, 100, 100);
    //    QList<Point> pts1 = getPoints(0, 10, 11, 11);
    //    QList<Point> feature1 = featurePointSelection(pts1, 11, 11);

    //    for (int i = 0; i < 100; i++)
    //    {
    //        cout << pts[i].x << "," << pts[i].y
    //             << "," << pts[i].z << endl;
    //    }
    //    cout << "===============================================" << endl;
    //    for (int i = 0; i < 100; i++)
    //    {
    //        cout << pts1[i].x << "," << pts1[i].y
    //             << "," << pts1[i].z << endl;
    //    }

    //    for (int i = 0; i < 10; i++)
    //    {
    //        for (int j = 0; j < 20; j++)
    //        {
    //            printf("%.1f ", pts[j*10+i].z);
    //        }
    //        printf("\n");
    //    }
    //    for (int i = 0; i < 10; i++)
    //    {
    //        for (int j = 0; j < 10; j++)
    //        {
    //            printf("%.1f ", pts1[j*10+i].z);
    //        }
    //        printf("\n");
    //    }

    //    cout << pts[0].x << "," << pts[0].y << endl;
    //    cout << pts[9900].x << "," << pts[9900].y << endl;
    //    cout << pts1[0].x << "," << pts1[0].y << endl;
    //    cout << pts1[9900].x << "," << pts1[9900].y << endl;

    //    QList<Point> pts2 = getPoints(100, 0, 101, 101);
    //    QList<Point> feature2 = featurePointSelection(pts2, 101, 101);
    //    QList<Point> pts3 = getPoints(100, 100, 101, 101);
    //    QList<Point> feature3 = featurePointSelection(pts3, 101, 101);
    //    int count = fearure.size();

    //    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    //    osg::Vec3dArray* points=new osg::Vec3dArray(count);
    //    for (int i=0;i<count;i++)
    //    {
    //        points->at(i)=osg::Vec3d(pts[i].x, pts[i].y, pts[i].z);
    //    }
    //    geom->setUseDisplayList(true);
    //    geom->setUseVertexBufferObjects(true);
    //    geom->setVertexArray(points);
    //    geom->addPrimitiveSet(new osg::DrawArrays(GL_TRIANGLE_FAN,0,count));
    //    cout << pts.size() << endl;
    //    cout << fearure.size() << endl;
    //    for (int i = 0; i < fearure1.size(); i++)
    //    {
    //        cout << fearure1[i].x << endl;
    //        fearure1[i].x -= 3000;
    //    }

    //    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
    //    osg::ref_ptr<osg::Group> scene=new osg::Group;
    ////    scene->addChild(geode.get());
    //    int interval = 100;
    //    for (int i = 0; i < 2; i++)
    //    {
    //        for (int j = 0; j < 2; j++)
    //        {
    //            QList<Point> pts = getPoints(i*interval, j*interval, interval+1, interval+1);
    ////            QList<Point> feature = featurePointSelection(pts, interval+1, interval+1);
    //            geode->addDrawable(createTriangulate(pts));
    ////            scene->addChild(createTriangulate(pts)->asTerrain());
    //        }
    //    }
    ////    geode->addDrawable(createTriangulate(feature));
    ////    geode->addDrawable(createTriangulate(feature1));
    ////    geode->addDrawable(createTriangulate(feature2));
    ////    geode->addDrawable(createTriangulate(feature3));
    ////    geode->addDrawable(createTriangulate(pts));
    ////    geode->addDrawable(createTriangulate(pts1));
    ////    geode->addDrawable(geom);

    //    osg::ref_ptr<osg::PolygonMode> pLolyMode= new osg::PolygonMode;
    //    osg::ref_ptr<osg::StateSet> stateset = new osg::StateSet;
    //    pLolyMode->setMode(osg::PolygonMode::FRONT_AND_BACK ,osg::PolygonMode::LINE);
    //    stateset->setAttribute(pLolyMode);
    //    geode->setStateSet(stateset);

    //    scene->addChild(geode.get());
    ////    scene->addChild(readFromDEM().get());

    //    osgUtil::Optimizer optimizer;
    //    optimizer.optimize(scene.get());

    //    osgViewer::Viewer viewer;
    //    viewer.setSceneData(scene);
    //    //        viewer.addEventHandler(new osgViewer::WindowSizeHandler());
    //    viewer.realize();
    //    return viewer.run();

    //    pts[100].print();
    //    pts[101].print();
    //    pts[102].print();
    //    pts[103].print();
    osg::ref_ptr<osg::Group> scene = new osg::Group;
    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
    osg::ref_ptr<osg::Node> node = new osg::Node();
//    node = osgDB::readNodeFile("terrain_no_normal.osg");
//    node = osgDB::readNodeFile("terrain.osg");
    node = osgDB::readNodeFile("/home/lzt/Downloads/osg-data/OpenSceneGraph-Data/cessna.osg");
//    osg::ref_ptr<osg::Image> image = osgDB::readImageFile("/home/lzt/material/img_map/gs_637_474_10.jpg");
//        osg::ref_ptr<osg::Image> image = osgDB::readImageFile("/home/lzt/Pictures/220956jq48448gomf05j4g.png");
        osg::ref_ptr<osg::Image> image = osgDB::readImageFile("/home/lzt/Pictures/IRTexture.bmp");
//        osg::ref_ptr<osg::Image> image = osgDB::readImageFile("/home/lzt/material/20150630152144255.png");

        osg::BoundingSphere bs;
        node->computeBound();
        bs = node->getBound();
        osg::ref_ptr<osg::Light> light = new osg::Light();
        light->setLightNum(0);
        //设置方向
        light->setDirection(osg::Vec3(0.0f, 0.0f, -1.0f));
        //设置位置
        light->setPosition(osg::Vec4(bs.center().x() + 600, bs.center().y(), bs.center().z()+ bs.radius(), 0.0f));
        //设置环境光的颜色
        light->setAmbient(osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));
        //设置散射光颜色
        light->setDiffuse(osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));

        //设置恒衰减指数
        light->setConstantAttenuation(1.0f);
        //设置线形衰减指数
        light->setLinearAttenuation(0.0f);
        //设置二次方衰减指数
//        light->setQuadraticAttenuation(0.0f);

        //创建光源
        osg::ref_ptr<osg::LightSource> lightSource = new osg::LightSource();
        lightSource->setLight(light.get());
//        scene->addChild(lightSource.get());


    if (image.get())
    {
        osg::ref_ptr<osg::Texture2D> texture=new osg::Texture2D();
        texture->setDataVariance(osg::Object::DYNAMIC);
        texture->setImage(image.get());
        texture->setUnRefImageDataAfterApply( true );
        //        texture->setWrap(osg::Texture2D::WrapParameter::WRAP_S,osg::Texture2D::WrapMode::CLAMP_TO_BORDER);
        //        texture->setWrap(osg::Texture2D::WrapParameter::WRAP_T,osg::Texture2D::WrapMode::CLAMP_TO_BORDER);
        //        texture->setWrap(osg::Texture2D::WrapParameter::WRAP_R,osg::Texture2D::WrapMode::CLAMP_TO_BORDER);

        //        texture->setFilter(osg::Texture2D::FilterParameter::MIN_FILTER,osg::Texture2D::FilterMode::LINEAR);
        //        texture->setFilter(osg::Texture2D::FilterParameter::MAG_FILTER,osg::Texture2D::FilterMode::LINEAR);

        node->getOrCreateStateSet()->setTextureAttributeAndModes(0, texture.get(), osg::StateAttribute::ON);
//        node->getOrCreateStateSet()->setTextureAttribute(0,texture,osg::StateAttribute::OVERRIDE);
//        node->getOrCreateStateSet()->setTextureMode(0,GL_TEXTURE_2D,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
//                node->getOrCreateStateSet()->setTextureMode(0,GL_TEXTURE_GEN_S,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
        //        node->getOrCreateStateSet()->setTextureMode(0,GL_TEXTURE_GEN_T,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
        ////                geometry->getOrCreateStateSet()->setAttribute(material,osg::StateAttribute::OVERRIDE);
    }
        node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);
        node->getOrCreateStateSet()->setMode(GL_LIGHT0, osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);
//    scene->addChild(node.get());
    scene->addChild(createTerrain());
    //    int xSize = 100, ySize = 100;

    //    for (int i = 0; i < 1; i++)
    //    {
    //        for (int j = 0; j < 1; j++)
    //        {
    //            vector<TRIANGLE_DESC> triangles;
    ////            QList<vec3> norms;
    //            osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
    //            triangles = getCVTriangles(i*xSize, j*ySize, xSize, ySize, true);
    ////            cout << "norms:" << norms.size() << endl;
    //            osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array();
    ////            osg::ref_ptr<osg::Vec3Array> n = new osg::Vec3Array();

    //            for (const auto triangle : triangles)
    //            {
    ////                osg::ref_ptr<osg::Vec3Array> n = new osg::Vec3Array();
    //                points->push_back(osg::Vec3(triangle.pt1.x, triangle.pt1.y, triangle.pt1.z));
    //                points->push_back(osg::Vec3(triangle.pt2.x, triangle.pt2.y, triangle.pt2.z));
    //                points->push_back(osg::Vec3(triangle.pt3.x, triangle.pt3.y, triangle.pt3.z));
    ////                n->push_back(osg::Vec3(triangle.norm1.x, triangle.norm1.y, triangle.norm1.z));
    ////                n->push_back(osg::Vec3(triangle.norm2.x, triangle.norm2.y, triangle.norm2.z));
    ////                n->push_back(osg::Vec3(triangle.norm3.x, triangle.norm3.y, triangle.norm3.z));
    //            }
    ////            for (int i = 0; i < triangles.size(); i++)
    ////            {
    ////                osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
    ////                osg::ref_ptr<osg::Vec3Array> smallpts;
    ////                osg::ref_ptr<osg::Vec3Array> smallnorm;
    ////                smallpts->push_back(points->at(i*3));
    ////                smallpts->push_back(points->at(i*3+1));
    ////                smallpts->push_back(points->at(i*3+2));
    ////                smallnorm->push_back(n->at(i));
    ////                geometry->setVertexArray(smallpts.get());
    ////                geometry->setNormalArray(smallnorm.get());
    ////                geometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
    ////                geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points->size()));
    ////                geode->addDrawable(geometry.get());
    ////            }
    //            geometry->setVertexArray(points.get());
    ////            geometry->setNormalArray(n.get());
    //            cout << "pt size" << points->size() << endl;
    ////            cout << "norm size" << n->size() << endl;
    ////            geometry->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
    //            geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points->size()));
    ////            osg::ref_ptr<osg::Vec2Array> vt = new osg::Vec2Array();

    ////            vt->push_back(osg::Vec2(0.0f, 0.0f));
    ////            vt->push_back(osg::Vec2(1.0f, 0.0f));
    ////            vt->push_back(osg::Vec2(1.0f, 1.0f));
    ////            vt->push_back(osg::Vec2(0.0f, 1.0f));

    ////            geometry->setTexCoordArray(0, vt.get());

    ////            osg::Texture2D* tex = new osg::Texture2D(image);

    ////            tex->setFilter(osg::Texture2D::MIN_FILTER,osg::Texture2D::LINEAR_MIPMAP_LINEAR);
    ////            tex->setFilter(osg::Texture2D::MAG_FILTER,osg::Texture2D::LINEAR);

    ////            tex->setWrap(osg::Texture::WRAP_S, osg::Texture::REPEAT);

    ////            tex->setWrap(osg::Texture::WRAP_T, osg::Texture::REPEAT);

    ////            geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0, tex);

    //            osgUtil::SmoothingVisitor smv;
    //            smv.smooth(*geometry);

    ////            osg::ref_ptr<osg::Vec2Array> texcoords = new osg::Vec2Array( 4 );
    ////            (*texcoords)[0].set( 0.0, 1.0 );
    ////            (*texcoords)[1].set( 0.0, 0.0 );
    ////            (*texcoords)[2].set( 1.0, 1.0 );
    ////            (*texcoords)[3].set( 1.0, 0.0 );
    ////            geometry->setTexCoordArray(0, texcoords.get());


    //            if (image.get())
    //            {
    //                osg::ref_ptr<osg::Texture2D> texture=new osg::Texture2D();
    //                texture->setDataVariance(osg::Object::DYNAMIC);
    //                texture->setImage(image.get());

    ////                texture->setWrap(osg::Texture2D::WrapParameter::WRAP_S,osg::Texture2D::WrapMode::CLAMP_TO_BORDER);
    ////                texture->setWrap(osg::Texture2D::WrapParameter::WRAP_T,osg::Texture2D::WrapMode::CLAMP_TO_BORDER);
    ////                texture->setWrap(osg::Texture2D::WrapParameter::WRAP_R,osg::Texture2D::WrapMode::CLAMP_TO_BORDER);

    ////                texture->setFilter(osg::Texture2D::FilterParameter::MIN_FILTER,osg::Texture2D::FilterMode::LINEAR);
    ////                texture->setFilter(osg::Texture2D::FilterParameter::MAG_FILTER,osg::Texture2D::FilterMode::LINEAR);

    //////                osg::material *material = new osg::material;
    //////                osg::stateset *stateset = new osg::stateset;

    ////                geode->getOrCreateStateSet()->setTextureAttribute(0,texture,osg::StateAttribute::OVERRIDE);
    ////                geode->getOrCreateStateSet()->setTextureMode(0,GL_TEXTURE_2D,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
    ////                geode->getOrCreateStateSet()->setTextureMode(0,GL_TEXTURE_GEN_S,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
    ////                geode->getOrCreateStateSet()->setTextureMode(0,GL_TEXTURE_GEN_T,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
    ////                geometry->getOrCreateStateSet()->setAttribute(material,osg::StateAttribute::OVERRIDE);


    ////                osg::ref_ptr<osg::TexGen> texgen = new osg::TexGen();
    ////                texgen->setMode(osg::TexGen::SPHERE_MAP);

    ////                osg::ref_ptr<osg::TexEnv> texenv = new osg::TexEnv();
    ////                texenv->setMode(osg::TexEnv::Mode::ADD);

    ////                geometry->getOrCreateStateSet()->setTextureAttributeAndModes(
    ////                            0, texture.get(), osg::StateAttribute::ON);
    ////                geometry->getOrCreateStateSet()->setTextureAttributeAndModes(
    ////                            0, texgen.get(), osg::StateAttribute::ON);
    ////                geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0, texenv.get());

    ////                //设置自动生成纹理坐标
    ////                osg::ref_ptr<osg::TexGen> texgen=new osg::TexGen();
    ////                texgen->setMode(osg::TexGen::NORMAL_MAP);

    ////                osg::ref_ptr<osg::TexEnv> texenv=new osg::TexEnv;
    ////                texenv->setMode(osg::TexEnv::DECAL);

    //                //启动单元一自动生成纹理坐标，并使用纹理
    ////                geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0,texture.get(),osg::StateAttribute::ON);
    ////                geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0,texgen.get(),osg::StateAttribute::ON);
    ////                geometry->getOrCreateStateSet()->setTextureAttribute(0,texenv.get());

    ////                texture->setFilter(osg::Texture2D::MIN_FILTER,osg::Texture2D::LINEAR_MIPMAP_LINEAR);
    ////                texture->setFilter(osg::Texture2D::MAG_FILTER,osg::Texture2D::LINEAR);

    ////                texture->setWrap(osg::Texture::WRAP_S, osg::Texture::REPEAT);
    ////                texture->setWrap(osg::Texture::WRAP_T, osg::Texture::REPEAT);
    ////                geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0, texture);
    //            }
    //            geode->addDrawable(geometry.get());

    //        }
    //    }


    //    geode->addDrawable(createTriangulate(pts));

    //    vector<TRIANGLE_DESC> triangles1;
    //    osg::ref_ptr<osg::Geometry> geometry1 = new osg::Geometry();
    //    triangles1 = getCVTriangles(0, 0, xSize, ySize, false);
    //    osg::ref_ptr<osg::Vec3Array> points1 = new osg::Vec3Array();

    //    for (const auto triangle : triangles1)
    //    {
    //        points1->push_back(osg::Vec3(triangle.pt1.x+3100, triangle.pt1.y, triangle.pt1.z));
    //        points1->push_back(osg::Vec3(triangle.pt2.x+3100, triangle.pt2.y, triangle.pt2.z));
    //        points1->push_back(osg::Vec3(triangle.pt3.x+3100, triangle.pt3.y, triangle.pt3.z));
    //    }
    //    geometry1->setVertexArray(points1.get());
    //    geometry1->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points1->size()));
    //    osgUtil::SmoothingVisitor smv1;
    //    smv1.smooth(*geometry1);

    //    geode->addDrawable(geometry1.get());


    //    osg::ref_ptr<osg::PolygonMode> polyMode= new osg::PolygonMode;
    //    osg::ref_ptr<osg::StateSet> stateset = scene->getOrCreateStateSet();
    //    polyMode->setMode(osg::PolygonMode::FRONT_AND_BACK ,osg::PolygonMode::LINE);
    //    stateset->setAttribute(polyMode);
    //    stateset->setMode(GL_LIGHTING, osg::StateAttribute::ON);
    //    stateset->setMode(GL_LIGHT0, osg::StateAttribute::ON);
    //    geode->setStateSet(stateset);

    //   /* osg::ref_ptr<osg::StateSet> */stateset = scene->getOrCreateStateSet();
    //    stateset->setMode(GL_LIGHTING, osg::StateAttribute::ON);
    //    stateset->setMode(GL_LIGHT0, osg::StateAttribute::ON);
    //    geode->setStateSet(stateset);

//    scene->addChild(geode.get());
    //    scene->addChild(readFromDEM().get());

    osgUtil::Optimizer optimizer;
    optimizer.optimize(scene.get());

//        osgDB::writeNodeFile(*(scene.get()), "terrain_no_normal.osg");

    osgViewer::Viewer viewer;
    viewer.setSceneData(scene);
    //        viewer.addEventHandler(new osgViewer::WindowSizeHandler());
    viewer.realize();
    return viewer.run();
}
