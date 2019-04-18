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
#include "demreader.h"
#include "skybox.h"

#define ROCK_TEXTURE 0
#define SAND_TEXTURE 1
#define TREE_TEXTURE 2
#define WATER_TEXTURE 3
#define NO_TEXTURE_GRAY 0
#define ROCK_TEXTURE_GRAY 0
#define SAND_TEXTURE_GRAY 86
#define TREE_TEXTURE_GRAY 171
#define WATER_TEXTURE_GRAY 255
#define ROCK_TEXTURE_B 0
#define ROCK_TEXTURE_G 0
#define ROCK_TEXTURE_R 255
#define SAND_TEXTURE_B 0
#define SAND_TEXTURE_G 255
#define SAND_TEXTURE_R 255
#define TREE_TEXTURE_B 0
#define TREE_TEXTURE_G 255
#define TREE_TEXTURE_R 0
#define WATER_TEXTURE_B 255
#define WATER_TEXTURE_G 0
#define WATER_TEXTURE_R 0

using namespace std;


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
    texture->setResizeNonPowerOfTwoHint(false);

    //创建状态集对象
    //    osg::ref_ptr<osg::StateSet> stateset = new osg::StateSet();
    //    stateset->setTextureAttributeAndModes(0, texture.get(), osg::StateAttribute::ON);

    return texture.release();
}

osg::ref_ptr<osg::Geode> createTerrain()
{
    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
//    osg::ref_ptr<osg::Image> image = new osg::Image; //= osgDB::readImageFile("/home/lzt/material/img_map/gs_637_474_10.jpg");
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array();
    osg::ref_ptr<osg::Vec2Array> texCoord = new osg::Vec2Array();
    DEMReader *demReader = new DEMReader("/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif");
    int startFromX = 3000, startFromY = 750;
    int xNum = 3, yNum = 3;
    int xSize = 100, ySize = 100;
    int zoomLevel = 15;
    bool infrared = true;
    double diffx = demReader->getDiffX();
    double diffy = demReader->getDiffY();
    double startLng = demReader->getOriginX() + startFromX*diffx;//start lng from the area you have chosen, has offset
    double startLat = demReader->getOriginY() + startFromY*diffy;
    double endLng = startLng + diffx*xSize*xNum;
    double endLat = startLat + diffy*ySize*yNum;
    cv::Point startTile;
    cv::Point endTile;
    demReader->getGoogleMapTile(startLng, startLat, &startTile, zoomLevel);
    demReader->getGoogleMapTile(endLng, endLat, &endTile, zoomLevel);
    cout << startTile.x << ", " << startTile.y << endl;
    cout << endTile.x << ", " << endTile.y << endl;
    int xPicNum = endTile.x-startTile.x;
    int yPicNum = endTile.y-startTile.y;
    cv::Mat bigImg(TILE_SIZE*(yPicNum+1), TILE_SIZE*(xPicNum+1), 16);

    for (int i = 0; i <= yPicNum; i++)
    {
        for (int j = 0; j <= xPicNum; j++)
        {
            QString imgFileName = "/home/lzt/material/img_map/";
            imgFileName.append(QString::number(zoomLevel));
            imgFileName.append("/gs_");
            imgFileName.append(QString::number(startTile.x+j));
            imgFileName.append("_");
            imgFileName.append(QString::number(startTile.y+i));
            imgFileName.append("_");
            imgFileName.append(QString::number(zoomLevel));
            imgFileName.append(".jpg");
            cv::Mat img = cv::imread(imgFileName.toStdString());
            cv::Mat tmp = bigImg.colRange(img.cols*j, img.cols*(j+1)).rowRange(img.rows*i, img.rows*(i+1));
            img.copyTo(tmp);
            cout << imgFileName.toStdString() << endl;
        }
    }
    if(infrared)
    {
        cv::Mat gray;
        cv::cvtColor(bigImg,gray,cv::COLOR_BGR2GRAY);
        cv::Mat materialMap = cv::imread("/home/lzt/material/img_map/tmp-IR.bmp", cv::IMREAD_GRAYSCALE);
        cout << materialMap.type() << endl;
        QVector<double> radiations;
        if ((materialMap.rows == gray.rows) && (materialMap.cols == gray.cols))
        {
            for (int i = 0; i < gray.rows; i++)
            {
                for (int j = 0; j < gray.cols; j++)
                {
                    int fixColor;
                    double fixRate;
                    switch ((int)materialMap.at<uchar>(i, j)) {
                    case ROCK_TEXTURE_GRAY:
                        fixColor = 50;
                        fixRate = 1.25;
                        break;
                    case SAND_TEXTURE_GRAY:
                        fixColor = 40;
                        fixRate = 1.2;
                        break;
                    case TREE_TEXTURE_GRAY:
                        fixColor = -70;
                        fixRate = 0.8;
                        break;
                    case WATER_TEXTURE_GRAY:
                        fixColor = -80;
                        fixRate = 0.75;
                        break;
                    default:
                        break;
                    }
                    int finalColor = (int)gray.at<uchar>(i, j) + fixColor;
//                    if (finalColor > 255)
//                        finalColor = 255;
//                    else if(finalColor < 0)
//                        finalColor = 0;
                    int finalRateColor = (int)gray.at<uchar>(i, j)*fixRate;
//                    if (finalRateColor == 0 || finalRateColor == 255)
//                        cout << i << ", " << j << endl;
//                    gray.at<uchar>(i, j) = finalRateColor;
                    radiations.push_back(finalColor);
                }
            }
            int minColor = 256, maxColor = -1;
            for (int i = 0; i < gray.rows; i++)
            {
                for (int j = 0; j < gray.cols; j++)
                {
//                    minColor = min(minColor, (int)gray.at<uchar>(i, j));
//                    maxColor = max(maxColor, (int)gray.at<uchar>(i, j));
                    minColor = min(minColor, (int)radiations.at(i*gray.cols+j));
                    maxColor = max(maxColor, (int)radiations.at(i*gray.cols+j));
                }
            }
            int intervals = maxColor - minColor;
            double gapColor = 255.0 / intervals;
            for (int i = 0; i < gray.rows; i++)
            {
                for (int j = 0; j < gray.cols; j++)
                {
                    gray.at<uchar>(i, j) = (radiations.at(i*gray.cols+j) - minColor)*gapColor;
                }
            }


        }
        cv::imwrite("tmp.jpg",gray);
    }else
    {
        cv::imwrite("tmp.jpg",bigImg);
    }


    int pixelStartFromX = startTile.x*TILE_SIZE;
    int pixelStartFromY = startTile.y*TILE_SIZE;

    for (int i = 0; i < yNum; i++)
    {
        for (int j = 0; j < xNum; j++)
        {
//            if (i==1&&j==1)
//                break;
            vector<TRIANGLE_DESC> triangles;
//            osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
            triangles = demReader->getCVTriangles(j*xSize + startFromX, i*ySize + startFromY, xSize+1, ySize+1, true);
//            osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array();
//            osg::ref_ptr<osg::Vec2Array> texCoord = new osg::Vec2Array();
            int height = demReader->getHeight();

            for (const auto triangle : triangles)
            {
                //insert the points reversely because that the y value has became an opposite of the original value
                points->push_back(osg::Vec3(triangle.pt3.x, (height-1)*30-triangle.pt3.y, triangle.pt3.z));
                points->push_back(osg::Vec3(triangle.pt2.x, (height-1)*30-triangle.pt2.y, triangle.pt2.z));
                points->push_back(osg::Vec3(triangle.pt1.x, (height-1)*30-triangle.pt1.y, triangle.pt1.z));
                //                double scale = triangle.pt3.x/(30*xSize);

                //calculate texture coordinates for Google map
                double lng, lat;
                int pixelX, pixelY;
                cv::Point pixelCoord;
                lng = startLng + ((triangle.pt3.x/30)-startFromX)*diffx;
                lat = startLat + ((triangle.pt3.y/30)-startFromY)*diffy;
                demReader->getGoogleMapPixel(lng, lat, &pixelCoord, zoomLevel);
                pixelX = pixelCoord.x - pixelStartFromX;
                pixelY = pixelCoord.y - pixelStartFromY;
                texCoord->push_back(osg::Vec2(pixelX/(bigImg.cols*1.0f), 1-pixelY/(bigImg.rows*1.0f)));

                lng = startLng + ((triangle.pt2.x/30)-startFromX)*diffx;
                lat = startLat + ((triangle.pt2.y/30)-startFromY)*diffy;
                demReader->getGoogleMapPixel(lng, lat, &pixelCoord, zoomLevel);
                pixelX = pixelCoord.x - pixelStartFromX;
                pixelY = pixelCoord.y - pixelStartFromY;
                texCoord->push_back(osg::Vec2(pixelX/(bigImg.cols*1.0f), 1-pixelY/(bigImg.rows*1.0f)));

                lng = startLng + ((triangle.pt1.x/30)-startFromX)*diffx;
                lat = startLat + ((triangle.pt1.y/30)-startFromY)*diffy;
                demReader->getGoogleMapPixel(lng, lat, &pixelCoord, zoomLevel);
                pixelX = pixelCoord.x - pixelStartFromX;
                pixelY = pixelCoord.y - pixelStartFromY;
                texCoord->push_back(osg::Vec2(pixelX/(bigImg.cols*1.0f), 1-pixelY/(bigImg.rows*1.0f)));
//                texCoord->push_back(osg::Vec2(triangle.pt3.x/(30*(xSize+1)), 1-triangle.pt3.y/(30*(ySize+1))));
//                texCoord->push_back(osg::Vec2(triangle.pt2.x/(30*(xSize+1)), 1-triangle.pt2.y/(30*(ySize+1))));
//                texCoord->push_back(osg::Vec2(triangle.pt1.x/(30*(xSize+1)), 1-triangle.pt1.y/(30*(ySize+1))));
            }
        }
    }
    osg::ref_ptr<osg::Image> image = osgDB::readImageFile("tmp.jpg");
//    image->setImage(bigImg.cols, bigImg.rows, 3, GL_BGR, GL_BGR, GL_UNSIGNED_BYTE, bigImg.data, osg::Image::NO_DELETE, 1);
    geometry->setVertexArray(points.get());
    geometry->setTexCoordArray(0,texCoord.get());
    cout << "pt size" << points->size() << endl;
    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, points->size()));
    geometry->getOrCreateStateSet()->setTextureAttributeAndModes(0, createTexture2D(image), osg::StateAttribute::ON);

    osgUtil::SmoothingVisitor smv;
    smv.smooth(*geometry);

    geode->addDrawable(geometry.get());
//    osg::ref_ptr<osg::PolygonMode> polyMode= new osg::PolygonMode;
//    polyMode->setMode(osg::PolygonMode::FRONT_AND_BACK ,osg::PolygonMode::LINE);
//    geode->getOrCreateStateSet()->setAttribute(polyMode);
    return geode.release();
}

int showTerrain()
{
    osg::ref_ptr<osg::Group> scene = new osg::Group;
//    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
//    osg::ref_ptr<osg::Node> node = osgDB::readNodeFile("terrain.osg");
//    osg::ref_ptr<osg::Node> model;
//    const float radius = 1.0f;
//    osg::Geode* geode = new osg::Geode;
//    geode->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(0.0f,0.0f,0.0f),radius)));
//    model = geode;
//    SkyBox *skybox = new SkyBox;
    scene->addChild(createTerrain());
//    scene->addChild(SkyBox::createSkyBox(true));
//    scene->addChild(model.get());
//    scene->addChild(node.get());

    osg::BoundingSphere bs;
    scene->getChild(0)->computeBound();
    bs = scene->getChild(0)->getBound();
    osg::ref_ptr<osg::Light> light = new osg::Light();
    light->setLightNum(0);
    //设置方向
    light->setDirection(osg::Vec3(0.0f, 0.0f, -1.0f));
    //设置位置
    light->setPosition(osg::Vec4(bs.center().x(), bs.center().y(), bs.center().z()+(bs.radius()), 0.0f));
    //设置环境光的颜色
    light->setAmbient(osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));
    //设置散射光颜色
    light->setDiffuse(osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));

    //设置恒衰减指数
    light->setConstantAttenuation(1.0f);
    //设置线形衰减指数
    light->setLinearAttenuation(0.0f);
//    //设置二次方衰减指数
//    //        light->setQuadraticAttenuation(0.0f);

//    //创建光源
    osg::ref_ptr<osg::LightSource> lightSource = new osg::LightSource();
    lightSource->setLight(light.get());
    scene->addChild(lightSource.get());

    osgUtil::Optimizer optimizer;
    optimizer.optimize(scene.get());

    osgViewer::Viewer viewer;
    viewer.setSceneData(scene);
    viewer.realize();
    return viewer.run();
}

int main()
{
    DEMReader *demReader = new DEMReader("/home/lzt/material/DEM/ASTGTM2_N12E044/ASTGTM2_N12E044_dem.tif");
    double diffx = demReader->getDiffX();
    double diffy = demReader->getDiffY();
    double startLng = demReader->getOriginX();
    double startLat = demReader->getOriginY();
    double endLng = startLng + diffx*3601;
    double endLat = startLat + diffy*3601;
    cv::Point startTile;
    cv::Point endTile;
    demReader->getGoogleMapTile(startLng, startLat, &startTile, 15);
    demReader->getGoogleMapTile(endLng, endLat, &endTile, 15);
    cout << startTile.x << ", " << startTile.y << endl;
    cout << endTile.x << ", " << endTile.y << endl;
    showTerrain();
//    #define CUBEMAP_FILENAME(face) "/home/lzt/material/sky_box/" #face ".jpg"
//    #define CUBEMAP_FILENAME_IR(face) "/home/lzt/material/sky_box/" #face "-ir.jpg"
//    cv::Mat matposx = cv::imread(CUBEMAP_FILENAME(posx));
//    cv::Mat matnegx = cv::imread(CUBEMAP_FILENAME(negx));
//    cv::Mat matposy = cv::imread(CUBEMAP_FILENAME(posy));
//    cv::Mat matnegy = cv::imread(CUBEMAP_FILENAME(negy));
//    cv::Mat matposz = cv::imread(CUBEMAP_FILENAME(posz));
//    cv::Mat matnegz = cv::imread(CUBEMAP_FILENAME(negz));
//    cv::Mat posxir, negxir, posyir, negyir, poszir, negzir;
//    cv::cvtColor(matposx, posxir, cv::COLOR_BGR2GRAY);
//    cv::imwrite(CUBEMAP_FILENAME_IR(posx), posxir);
//    cv::cvtColor(matnegx, negxir, cv::COLOR_BGR2GRAY);
//    cv::imwrite(CUBEMAP_FILENAME_IR(negx), negxir);
//    cv::cvtColor(matposy, posyir, cv::COLOR_BGR2GRAY);
//    cv::imwrite(CUBEMAP_FILENAME_IR(posy), posyir);
//    cv::cvtColor(matnegy, negyir, cv::COLOR_BGR2GRAY);
//    cv::imwrite(CUBEMAP_FILENAME_IR(negy), negyir);
//    cv::cvtColor(matposz, poszir, cv::COLOR_BGR2GRAY);
//    cv::imwrite(CUBEMAP_FILENAME_IR(posz), poszir);
//    cv::cvtColor(matnegz, negzir, cv::COLOR_BGR2GRAY);
//    cv::imwrite(CUBEMAP_FILENAME_IR(negz), negzir);

}
