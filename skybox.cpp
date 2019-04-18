#include "skybox.h"
#include <osgDB/ReadFile>
#include <osg/TexGen>
#include <osg/TexEnv>
#include <osg/Depth>
#include <osg/ShapeDrawable>

SkyBox::SkyBox()
{

}

osg::TextureCubeMap* SkyBox::readCubeMap(bool infrared)
{
    osg::TextureCubeMap* cubemap = new osg::TextureCubeMap;
    //#define CUBEMAP_FILENAME(face) "nvlobby_" #face ".png"
    //#define CUBEMAP_FILENAME(face) "Cubemap_axis/" #face ".png"
//    #define CUBEMAP_FILENAME(face) "/home/lzt/Downloads/osg-data/OpenSceneGraph-Data/Cubemap_axis/" #face ".png"
#define CUBEMAP_FILENAME(face) "/home/lzt/material/sky_box/" #face ".jpg"
#define CUBEMAP_FILENAME_IR(face) "/home/lzt/material/sky_box/" #face "-ir.jpg"

    osg::ref_ptr<osg::Image>imagePosX;
    osg::ref_ptr<osg::Image>imageNegX;
    osg::ref_ptr<osg::Image>imagePosY;
    osg::ref_ptr<osg::Image>imageNegY;
    osg::ref_ptr<osg::Image>imagePosZ;
    osg::ref_ptr<osg::Image>imageNegZ;

    if (infrared)
    {
        imagePosX = osgDB::readRefImageFile(CUBEMAP_FILENAME_IR(posx));
        imageNegX = osgDB::readRefImageFile(CUBEMAP_FILENAME_IR(negx));
        imagePosY = osgDB::readRefImageFile(CUBEMAP_FILENAME_IR(posy));
        imageNegY = osgDB::readRefImageFile(CUBEMAP_FILENAME_IR(negy));
        imagePosZ = osgDB::readRefImageFile(CUBEMAP_FILENAME_IR(posz));
        imageNegZ = osgDB::readRefImageFile(CUBEMAP_FILENAME_IR(negz));
    }else
    {
        imagePosX = osgDB::readRefImageFile(CUBEMAP_FILENAME(posx));
        imageNegX = osgDB::readRefImageFile(CUBEMAP_FILENAME(negx));
        imagePosY = osgDB::readRefImageFile(CUBEMAP_FILENAME(posy));
        imageNegY = osgDB::readRefImageFile(CUBEMAP_FILENAME(negy));
        imagePosZ = osgDB::readRefImageFile(CUBEMAP_FILENAME(posz));
        imageNegZ = osgDB::readRefImageFile(CUBEMAP_FILENAME(negz));
    }


    if (imagePosX && imageNegX && imagePosY && imageNegY && imagePosZ && imageNegZ)
    {
        cubemap->setImage(osg::TextureCubeMap::POSITIVE_X, imagePosX);
        cubemap->setImage(osg::TextureCubeMap::NEGATIVE_X, imageNegX);
        cubemap->setImage(osg::TextureCubeMap::POSITIVE_Y, imagePosY);
        cubemap->setImage(osg::TextureCubeMap::NEGATIVE_Y, imageNegY);
        cubemap->setImage(osg::TextureCubeMap::POSITIVE_Z, imagePosZ);
        cubemap->setImage(osg::TextureCubeMap::NEGATIVE_Z, imageNegZ);

        cubemap->setWrap(osg::Texture::WRAP_S, osg::Texture::CLAMP_TO_EDGE);
        cubemap->setWrap(osg::Texture::WRAP_T, osg::Texture::CLAMP_TO_EDGE);
        cubemap->setWrap(osg::Texture::WRAP_R, osg::Texture::CLAMP_TO_EDGE);

        cubemap->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);
        cubemap->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
    }
    else
    {
        std::cout << "no img!" << std::endl;
        std::cout << CUBEMAP_FILENAME(posx) << std::endl;
    }

    return cubemap;
}

osg::Node *SkyBox::createSkyBox(bool infrared)
{

    osg::StateSet* stateset = new osg::StateSet();

    osg::TexEnv* te = new osg::TexEnv;
    te->setMode(osg::TexEnv::REPLACE);
    stateset->setTextureAttributeAndModes(0, te, osg::StateAttribute::ON);

    osg::TexGen *tg = new osg::TexGen;
    tg->setMode(osg::TexGen::NORMAL_MAP);
    stateset->setTextureAttributeAndModes(0, tg, osg::StateAttribute::ON);

    osg::TexMat *tm = new osg::TexMat;
    stateset->setTextureAttribute(0, tm);

    osg::TextureCubeMap* skymap = readCubeMap(infrared);
    stateset->setTextureAttributeAndModes(0, skymap, osg::StateAttribute::ON);

    stateset->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
    stateset->setMode( GL_CULL_FACE, osg::StateAttribute::OFF );

    // clear the depth to the far plane.
    osg::Depth* depth = new osg::Depth;
    depth->setFunction(osg::Depth::ALWAYS);
    depth->setRange(1.0,1.0);
    stateset->setAttributeAndModes(depth, osg::StateAttribute::ON );

    stateset->setRenderBinDetails(-1,"RenderBin");

//    int startFromX = 3000, startFromY = 750;
//    int xNum = 3, yNum = 3;
//    int xSize = 100, ySize = 100;
//    double centerx = (startFromX + (xSize+xNum)/2.0)*30;
//    double centery = (startFromY + (ySize+yNum)/2.0)*30;

    osg::Drawable* drawable = new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(0.0f,0.0f,0.0f),1000));

    osg::Geode* geode = new osg::Geode;
    geode->setCullingActive(false);
    geode->setStateSet( stateset );
    geode->addDrawable(drawable);


    osg::Transform* transform = new MoveEarthySkyWithEyePointTransform;
    transform->setCullingActive(false);
    transform->addChild(geode);

    osg::ClearNode* clearNode = new osg::ClearNode;
//  clearNode->setRequiresClear(false);
    clearNode->setCullCallback(new TexMatCallback(*tm));
    clearNode->addChild(transform);

    return clearNode;
}
