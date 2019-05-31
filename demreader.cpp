#include "demreader.h"
using namespace std;

DEMReader::DEMReader(string filename, double theta_threshold, double phi_threshold)
{
//    const char *pszFilename="/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif";
    this->theta_threshold = theta_threshold;
    this->phi_threshold = phi_threshold;
    GDALAllRegister();
    poDataset = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
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
    this->width = poDataset->GetRasterXSize();
    this->height = poDataset->GetRasterYSize();

//    GDALRasterBand *poBand;
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
}

DEMReader::~DEMReader()
{
    GDALClose((GDALDatasetH)poDataset);
}

cv::Point2d DEMReader::project(double lat, double lng) {
    double siny = sin(lat * M_PI / 180);

    // Truncating to 0.9999 effectively limits latitude to 89.189. This is
    // about a third of a tile past the edge of the world tile.
    siny = min(max(siny, -0.9999), 0.9999);
    double x = TILE_SIZE * (0.5 + lng / 360);
    double y = TILE_SIZE * (0.5 - log((1 + siny) / (1 - siny)) / (4 * M_PI));

    return cv::Point2d(x, y);
}

QList<Point> DEMReader::getPoints(int startX, int startY, int sizeX, int sizeY)
{
    //    const char *pszFilename="/home/lzt/apps/tin-terrain/3rdparty/craterlake/dems_10m.dem";
    const char *pszFilename="/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif";
    //"/home/lzt/material/DEM/n39_e126_1arc_v3.tif";
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
    width = poDataset->GetRasterXSize();
    height = poDataset->GetRasterYSize();

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

osg::ref_ptr<osg::Geode> DEMReader::getTerrain(int startX, int startY, int sizeX, int sizeY)
{
    this->startX = startX;
    this->startY = startY;
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    QList<Point> pts;
    QVector<double> elev;

    float *pafScanline;
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
            x = (j+startX)*30;
            y = (i+startY)*30;
            pts.push_back(Point(x, y, elevation));
            elev.push_back(elevation);
            pafScanline++;
        }
    }
        QList<Point> feature = featurePointSelection(pts, sizeX, sizeY/*, resNorms*/);
        osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array();
        osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array();
        for (const auto point : feature)
        {
            coords->push_back(osg::Vec3(point.x, point.y, point.z));
        }

        osg::ref_ptr<osgUtil::DelaunayTriangulator> dt = new osgUtil::DelaunayTriangulator();
        dt->setInputPointArray(coords);//赋给它三维点集数组
        dt->setOutputNormalArray(normals);//输出法向量
        dt->triangulate();
        osg::ref_ptr<deprecated_osg::Geometry> geometry = new deprecated_osg::Geometry();
        geometry->setVertexArray(coords.get());
        geometry->addPrimitiveSet(dt->getTriangles());
        geometry->setNormalArray(normals.get());
//        Geometry::fixDeprecatedData();
        geometry->setNormalBinding(deprecated_osg::Geometry::AttributeBinding::BIND_PER_PRIMITIVE);
        osg::ref_ptr<osg::Geode> geode = new osg::Geode();
        geode->addDrawable(geometry.get());

        return geode;
}

vector<TRIANGLE_DESC> DEMReader::getCVTriangles(int startX, int startY, int sizeX, int sizeY, bool needFeature)
{
    this->startX = startX;
    this->startY = startY;
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    QList<Point> pts;
    QVector<double> elev;

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
//            y = (originY-nYSize-i-1)*30;
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
//    double lng = adfGeoTransform[0]+(startX)*adfGeoTransform[1]+(startY)*adfGeoTransform[2];
//    double lat = adfGeoTransform[3]+(startX)*adfGeoTransform[4]+(startY)*adfGeoTransform[5];
//    cv::Point2d worldCoord = project(lat, lng);
//    cout << "world:" << worldCoord.x << ", " << worldCoord.y << endl;
//    int scale = 1 << 12;
//    cv::Point pixelCoord = cv::Point(floor(worldCoord.x * scale), floor(worldCoord.y * scale));
//    cout << "pixel:" << pixelCoord.x << ", " << pixelCoord.y << endl;
//    cv::Point tileCoord = cv::Point(floor(worldCoord.x * scale/TILE_SIZE), floor(worldCoord.y * scale/TILE_SIZE));
//    cout << "tile:" << tileCoord.x << ", " << tileCoord.y << endl;
//    cout << "pixel in image:" << pixelCoord.x-tileCoord.x*TILE_SIZE << ", " << pixelCoord.y-tileCoord.y*TILE_SIZE << endl;

//    lng = adfGeoTransform[0]+(startX+sizeX)*adfGeoTransform[1]+(startY+sizeY)*adfGeoTransform[2];
//    lat = adfGeoTransform[3]+(startX+sizeX)*adfGeoTransform[4]+(startY+sizeY)*adfGeoTransform[5];
//    lng = adfGeoTransform[0]+(3601)*adfGeoTransform[1]+(3601)*adfGeoTransform[2];
//    lat = adfGeoTransform[3]+(3601)*adfGeoTransform[4]+(3601)*adfGeoTransform[5];
//    worldCoord = project(lat, lng);
//    cout << "world:" << worldCoord.x << ", " << worldCoord.y << endl;
//    pixelCoord = cv::Point(floor(worldCoord.x * scale), floor(worldCoord.y * scale));
//    cout << "pixel:" << pixelCoord.x << ", " << pixelCoord.y << endl;
//    tileCoord = cv::Point(floor(worldCoord.x * scale/TILE_SIZE), floor(worldCoord.y * scale/TILE_SIZE));
//    cout << "tile:" << tileCoord.x << ", " << tileCoord.y << endl;
//    cout << "pixel in image:" << pixelCoord.x-tileCoord.x*TILE_SIZE << ", " << pixelCoord.y-tileCoord.y*TILE_SIZE << endl;

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

vector<TRIANGLE_DESC> DEMReader::getCVTrianglesNFeature(int startX, int startY, int sizeX, int sizeY, bool needFeature, QList<Point> *feature)
{
    this->startX = startX;
    this->startY = startY;
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    QList<Point> pts;
    QVector<double> elev;

    float *pafScanline;
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
            x = (j+startX)*30;
            y = (i+startY)*30;
            pts.push_back(Point(x, y, elevation));
            elev.push_back(elevation);
            pafScanline++;
        }
    }

    *feature = featurePointSelection(pts, sizeX, sizeY/*, resNorms*/);

    //opencv delaunay
    const int width = sizeX*32;
    const int height = sizeY*32;
    double offX, offY;
    offX = startX*30;
    offY = startY*30;
    vector<cv::Point3d> testPoints;
    if (needFeature)
    {
        for (int i = 0; i < feature->size(); i++)
        {
            testPoints.push_back(cv::Point3d(feature->at(i).x, feature->at(i).y, feature->at(i).z));
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

QList<Point> DEMReader::featurePointSelection(QList<Point> pts, int row, int col)
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
//                if (j == 0 || j == col-1)
//                {
//                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                }
//                else
//                {
//                    double diff1, diff2, diff3;
//                    //same or different direction of slope means different gradient
//                    diff1 = pts[i*col+j].z - pts[i*col+j-1].z;
//                    diff2 = pts[i*col+j+1].z - pts[i*col+j].z;
//                    diff3 = fabs(diff1-diff2);
//                    //not the corner point
//                    if (diff3 > SIDE_THRESHOLD)
//                    {
//                        respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                    }else if (i == 0 && pts[i*col+j].z == 0 && (pts[i*col+j-1].z != 0 || pts[i*col+j+1].z != 0
//                                                                || pts[(i+1)*col+j].z != 0 || pts[(i+1)*col+j-1].z != 0
//                                                                || pts[(i+1)*col+j+1].z != 0))
//                    {
//                        respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                    }else if (i == row-1 && pts[i*col+j].z == 0 && (pts[i*col+j-1].z != 0 || pts[i*col+j+1].z != 0
//                                                                    || pts[(i-1)*col+j].z != 0 || pts[(i-1)*col+j-1].z != 0
//                                                                    || pts[(i-1)*col+j+1].z != 0))
//                    {
//                        respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                    }

//                }

                respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
            }
            //for the first and last column(without cornor points)
            else if (j == 0 || j == col-1)
            {
//                double diff1, diff2, diff3;
//                //same or different direction of slope means different gradient
//                diff1 = pts[i*col+j].z - pts[(i-1)*col+j].z;
//                diff2 = pts[(i+1)*col+j].z - pts[i*col+j].z;
//                diff3 = fabs(diff1-diff2);
//                if (diff3 > SIDE_THRESHOLD)
//                {
//                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                }else if (j == 0 && pts[i*col+j].z == 0 && (pts[(i-1)*col+j].z != 0 || pts[(i+1)*col+j].z != 0
//                                                            || pts[i*col+j+1].z != 0 || pts[(i-1)*col+j+1].z != 0
//                                                            || pts[(i+1)*col+j+1].z != 0))
//                {
//                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                }else if (j == col-1 && pts[i*col+j].z == 0 && (pts[(i-1)*col+j].z != 0 || pts[(i+1)*col+j].z != 0
//                                                                || pts[i*col+j-1].z != 0 || pts[(i+1)*col+j-1].z != 0
//                                                                || pts[(i-1)*col+j-1].z != 0))
//                {
//                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
//                }
                respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
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

//                if (i == 4 && j == 134)
//                {
//                    cout << "here" << endl;
//                }
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
                    if (diffThetas[k] > theta_threshold)
                        count++;
                    if (diffPhis[k] > phi_threshold)
                        count++;
                }
                if (count > MAX_COUNT_NUM)
                {
                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
                }/*else if(((pts[i*col+j].z == pts[i*col+j-1].z)&&(pts[i*col+j].z != pts[i*col+j+1].z)) ||
                         (pts[i*col+j].z == pts[i*col+j+1].z)&&(pts[i*col+j].z != pts[i*col+j-1].z))
                {
                    respts.push_back(Point(pts[i*col+j].x, pts[i*col+j].y, pts[i*col+j].z));
                }*/else if (pts[i*col+j].z == 0 && (pts[i*col+j-1].z != 0 || pts[i*col+j+1].z != 0
                                                    || pts[(i-1)*col+j].z != 0 || pts[(i+1)*col+j].z != 0
                                                    || pts[(i-1)*col+j-1].z != 0 || pts[(i-1)*col+j+1].z != 0
                                                    || pts[(i+1)*col+j-1].z != 0 || pts[(i+1)*col+j+1].z != 0))
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

vector<TRIANGLE_DESC> DEMReader::delaunayAlgorithm(const cv::Rect boundRc,const vector<cv::Point3d>& points, QVector<double> elev, int col, int xOff, int yOff)
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

void DEMReader::calTextureRange()
{

}

int DEMReader::getHeight()
{
    return height;
}

int DEMReader::getWidth()
{
    return width;
}

double DEMReader::getOriginX()
{
    return originX;
}

double DEMReader::getOriginY()
{
    return originY;
}

double DEMReader::getDiffX()
{
    return diffX;
}
double DEMReader::getDiffY()
{
    return diffY;
}

void DEMReader::getGoogleMapPixel(double lng, double lat, cv::Point *pixelCoord, int zoomLevel)
{
    cv::Point2d worldCoord = project(lat, lng);
//    cout << "world:" << worldCoord.x << ", " << worldCoord.y << endl;
    int scale = 1 << zoomLevel;
    *pixelCoord = cv::Point(floor(worldCoord.x * scale), floor(worldCoord.y * scale));
//    cout << "pixel:" << pixelCoord.x << ", " << pixelCoord.y << endl;
//    *tileCoord = cv::Point(floor(worldCoord.x * scale/TILE_SIZE), floor(worldCoord.y * scale/TILE_SIZE));
//    cout << "tile:" << tileCoord.x << ", " << tileCoord.y << endl;
//    cout << "pixel in image:" << pixelCoord.x-tileCoord.x*TILE_SIZE << ", " << pixelCoord.y-tileCoord.y*TILE_SIZE << endl;
}

void DEMReader::getGoogleMapTile(double lng, double lat, cv::Point *tileCoord, int zoomLevel)
{
    cv::Point2d worldCoord = project(lat, lng);
//    cout << "world:" << worldCoord.x << ", " << worldCoord.y << endl;
    int scale = 1 << zoomLevel;
//    cv::Point *pixelCoord = cv::Point(floor(worldCoord.x * scale), floor(worldCoord.y * scale));
//    cout << "pixel:" << pixelCoord.x << ", " << pixelCoord.y << endl;
    *tileCoord = cv::Point(floor(worldCoord.x * scale/TILE_SIZE), floor(worldCoord.y * scale/TILE_SIZE));
//    cout << "tile:" << tileCoord.x << ", " << tileCoord.y << endl;
//    cout << "pixel in image:" << pixelCoord.x-tileCoord.x*TILE_SIZE << ", " << pixelCoord.y-tileCoord.y*TILE_SIZE << endl;
}

//osg::Drawable *createTriangulate(QList<Point> pts)
//{
//    //    osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array();
//    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array();
//    //    osg::ref_ptr<osg::Vec3Array> norms = new osg::Vec3Array();
//    //    float vertex[][3] = { -5.0f,-5.0f, 0.4f,
//    //                          1.0f, -5.6f, 0.0f,
//    //                          5.0f, -4.0f, -0.5f,
//    //                          -6.2f, 0.0f, 4.2f,
//    //                          -1.0f,-0.5f, 4.8f,
//    //                          4.3f, 1.0f, 3.0f,
//    //                          -4.8f, 5.4f, 0.3f,
//    //                          0.6f, 5.1f,-0.8f,
//    //                          5.2f, 4.5f, 0.1f };
//    //    unsigned int n = sizeof(vertex) / sizeof(float[3]);
//    for (unsigned int i = 0; i < pts.size(); i++)
//    {
//        //        coords->push_back(osg::Vec3(pts[i].x, pts[i].y, 0));
//        points->push_back(osg::Vec3(pts[i].x, pts[i].y, pts[i].z));
//    }
//    cout << pts.size() << "==" << points->size() << endl;
//    cout << "converting..." << endl;
//    osg::ref_ptr<osgUtil::DelaunayTriangulator> dt = new osgUtil::DelaunayTriangulator();
//    dt->setInputPointArray(points);
//    //    dt->setOutputNormalArray(norms);
//    dt->triangulate();
//    cout << "converting complete" << endl;
//    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
//    geometry->setVertexArray(points.get());
//    //    geometry->setNormalArray(norms.get());
//    //    geometry->setNormalBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);

//    geometry->addPrimitiveSet(dt->getTriangles());
//    //    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points->size()));
//    return geometry.release();

//}
