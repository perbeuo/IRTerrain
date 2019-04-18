#ifndef MYHANDLER_H
#define MYHANDLER_H

#include <QMap>
#include <iomanip>
#include <iostream>
#include "stdlib.h"
#include <osg/MatrixTransform>
#include <osgGA/GUIActionAdapter>
#include <osgGA/GUIEventHandler>
#include <osgGA/TrackballManipulator>
#include <osgSim/DOFTransform>
#include <osgText/Text>
#include <osgViewer/Viewer>
#include <osgUtil/SceneView>
#include <sstream>

// similar to TrackballManipulator(keydown and frame handler
// change)
class MyHandler : public osgGA::TrackballManipulator
{
public:
    /**
     * @brief MyHandler
     */
    MyHandler();
    ~MyHandler() {}

    /**
     * @brief handleKeyDown
     * @param ea
     * @param aa
     * @return
     */
    virtual bool handle(
            const osgGA::GUIEventAdapter& ea,
            osgGA::GUIActionAdapter& aa);

//    /**
//     * @brief setLabel
//     * @param name
//     */
//    void setLabel(const std::string& name);

    /**
     * @brief performMovementLeftKey
     * @param eventTimeDelta
     * @param dx
     * @param dy
     * @return
     */
    bool performMovementLeftKey(
            double eventTimeDelta, double dx, double dy);

    void rotateCameraLR(double theta, osgViewer::Viewer *viewer);

//    void rotateCameraUP(double val);

    void rotateCameraUD(double theta, osgViewer::Viewer *viewer);

//    virtual osg::Matrixd getInverseMatrix() const;

    osg::Vec3d world2Screen(osg::Vec3 worldPoint, osgViewer::Viewer *viewer);//世界到屏幕

    osg::Vec3d screen2World(osg::Vec3 screenPoint, osgViewer::Viewer* viewer);

    osg::Vec3d world2Camera(osg::Vec3 worldPoint, osgViewer::Viewer* viewer);//世界转相机

    osg::Vec3d camera2World(osg::Vec3 cameraPoint, osgViewer::Viewer* viewer);//相机转世界

protected:
    // HUD数据栏文本节点
    osg::ref_ptr<osgText::Text> text;

    osg::Vec3	_eye;				//视点位置
    osg::Vec3   _direction;         //视点方向
    osg::Vec3	_up;
};

#endif  // MYHANDLER_H
