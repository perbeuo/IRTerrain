#include "myhandler.h"
#include <QDebug>
using std::endl;

// const
//static const double ANGLE = 180 / 3.1415926;
//static const double TEMP = 180 / 20037508.34;
//extern std::map<unsigned long, MpsSimu::MpsIDLUAV> uavs;

MyHandler::MyHandler()
{
    computeHomePosition();
    _eye = _homeEye;
    _direction = _homeCenter - _homeEye;
    _up = _homeUp;
}

//osg::Matrixd MyHandler::getInverseMatrix() const
//{
//    osg::Matrix mat;
//    mat.makeLookAt(_eye, _eye + _direction, _up);

//    return mat;
//}


bool MyHandler::handle(
        const osgGA::GUIEventAdapter& eventAdapter,
        osgGA::GUIActionAdapter& actionAdapter)
{
    osgViewer::Viewer* viewer =
            dynamic_cast<osgViewer::Viewer*>(
                &actionAdapter);
    if (viewer == nullptr) return false;
    switch (eventAdapter.getEventType()) {
    case osgGA::GUIEventAdapter::FRAME:
        StandardManipulator::handleFrame(
                    eventAdapter, actionAdapter);
        break;

    case osgGA::GUIEventAdapter::KEYDOWN:

        //空格键使得相机回到home处
        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_Space)
        {
            flushMouseEventStack();
            _thrown = false;
            home(eventAdapter, actionAdapter);
            return true;
        }//COOR_CONVERT_FACTOR

        //键盘A键
        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_A)
        {
            flushMouseEventStack();
            _thrown = false;
            panModel(-5.0f, .0f, .0f);
            actionAdapter.requestRedraw();
        }

        //键盘D键
        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_D)
        {
            flushMouseEventStack();
            _thrown = false;
            panModel(5.0f, .0f, .0f);
            actionAdapter.requestRedraw();
        }

        //键盘W键
        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_W)
        {
            flushMouseEventStack();
            _thrown = false;
            panModel(.0f, .0f, -5.0f);
            actionAdapter.requestRedraw();
        }

        //键盘S键
        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_S)
        {
            flushMouseEventStack();
            _thrown = false;
            panModel(.0f, .0f, 5.0f);
            actionAdapter.requestRedraw();
        }

        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_Up)
        {
            flushMouseEventStack();
            _thrown = false;
            rotateCameraUD(-.25f, viewer);
            actionAdapter.requestRedraw();
        }

        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_Down)
        {
            flushMouseEventStack();
            _thrown = false;
            rotateCameraUD(.25f, viewer);
            actionAdapter.requestRedraw();
        }

        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_Left)
        {
            flushMouseEventStack();
            _thrown = false;
            rotateCameraLR(-.5f, viewer);
            actionAdapter.requestRedraw();
        }

        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_Right)
        {
            flushMouseEventStack();
            _thrown = false;
            rotateCameraLR(.5f, viewer);
            actionAdapter.requestRedraw();
        }

        if (eventAdapter.getKey() ==
                osgGA::GUIEventAdapter::KEY_P)
        {
            flushMouseEventStack();
            _thrown = false;

            osg::Vec3d eye,center,up;
            getTransformation(eye, center,up);

            QString msg0 =
                    "eye x: " + QString::number(eye.x())
                    + " y: " + QString::number(eye.y())
                    + " z: " + QString::number(eye.z());
            qDebug() << msg0 << endl;
            QString msg1 =
                    "center x: " + QString::number(center.x())
                    + " y: " + QString::number(center.y())
                    + " z: " + QString::number(center.z());
            qDebug() << msg1 << endl;
            QString msg2 =
                    "up x: " + QString::number(up.x())
                    + " y: " + QString::number(up.y())
                    + " z: " + QString::number(up.z());
            qDebug() << msg2 << endl;
            return true;
        }
        break;

    case osgGA::GUIEventAdapter::DRAG:
        StandardManipulator::handleMouseDrag(
                    eventAdapter, actionAdapter);
        break;
    case osgGA::GUIEventAdapter::PUSH:
        StandardManipulator::handleMousePush(
                    eventAdapter, actionAdapter);
        break;
    case osgGA::GUIEventAdapter::RELEASE:
        StandardManipulator::handleMouseRelease(
                    eventAdapter, actionAdapter);
        break;
    case osgGA::GUIEventAdapter::RESIZE:
        StandardManipulator::handleResize(
                    eventAdapter, actionAdapter);
        break;
    case osgGA::GUIEventAdapter::SCROLL:
        if (eventAdapter.getScrollingMotion() == osgGA::GUIEventAdapter::SCROLL_UP)
        {
            panModel(.0f, .0f, 200.0f);
        }
        else if (eventAdapter.getScrollingMotion() == osgGA::GUIEventAdapter::SCROLL_DOWN)
        {
            panModel(.0f, .0f, -200.0f);
        }

        break;
    default:
        break;
    }

    return false;
}

//void MyHandler::setLabel(const std::string& name)
//{
//    if (text.get()) text->setText(name);
//}

//键盘ADWS键控制视角的yaw、pitch角
bool MyHandler::performMovementLeftKey(
        double eventTimeDelta, double dx, double dy)
{
    // rotate camera
    if (getVerticalAxisFixed())
        rotateWithFixedVertical(dx, dy);
    else
    {
        rotateTrackball(
                    0.f, 0.f, dx, dy,
                    getThrowScale(eventTimeDelta));
    }

    return true;
}

void MyHandler::rotateCameraLR(double theta, osgViewer::Viewer *viewer)
{
    osg::Vec3d eye,center,up;
    getTransformation(eye, center,up);
    osg::Vec3d object, cam, cam1;
    object.x() = center.x();
    object.y() = center.y();
    object.z() = center.z();
    cam = world2Camera(object, viewer);

    //right-left
    double r = sqrt(pow(cam.x(),2)
                    +pow(cam.z(),2));
    double base_angle =
            atan2(cam.z(),
                  cam.x())*180/M_PI;
    double newX = r *
            cos((base_angle + theta) * M_PI / 180);
    double newZ = r *
            sin((base_angle + theta) * M_PI / 180);
    cam1.x() = newX;//forward-backward
    cam1.y() = cam.y();//left-right
    cam1.z() = newZ;//up-down
    setHomePosition(eye, camera2World(cam1, viewer), up);
    home(0);
}

//void MyHandler::rotateCameraUP(double val)
//{
//    osg::Vec3d eye,center,up,up1;
//    getTransformation(eye, center,up);

//    up1.x() = up.x()+ val;
//    up1.y() = up.y()+ val;
//    up1.z() = up.z();
//    setHomePosition(eye, center, up1);
//    home(0);
//}

void MyHandler::rotateCameraUD(double theta, osgViewer::Viewer* viewer){
    osg::Vec3d eye,center,up;
    getTransformation(eye, center,up);
    osg::Vec3d object, cam, cam1;
    object.x() = center.x();
    object.y() = center.y();
    object.z() = center.z();
    cam = world2Camera(object, viewer);

    double r = sqrt(pow(cam.z(),2)
                    +pow(cam.y(),2));
    double base_angle =
            atan2(cam.y(),
                  cam.z())*180/M_PI;
    double newZ = r *
            cos((base_angle + theta) * M_PI / 180);
    double newY = r *
            sin((base_angle + theta) * M_PI / 180);
    cam1.x() = cam.x();//forward-backward
    cam1.y() = newY;//left-right
    cam1.z() = newZ;//up-down
    setHomePosition(eye, camera2World(cam1, viewer), up);
    home(0);
}

osg::Vec3d MyHandler::world2Screen(osg::Vec3 worldPoint, osgViewer::Viewer* viewer)//世界到屏幕
{
    osg::Vec3d vec3;
    osg::ref_ptr<osg::Camera> camera = viewer->getCamera();
    osg::Matrix mVPW = camera->getViewMatrix() * camera->getProjectionMatrix() * camera->getViewport()->computeWindowMatrix();
    vec3 = worldPoint * mVPW;
    return vec3;
}

osg::Vec3d MyHandler::screen2World(osg::Vec3 screenPoint, osgViewer::Viewer* viewer)//将屏幕坐标转换到世界坐标
{
    osg::Vec3d vec3;
    osg::ref_ptr<osg::Camera> camera = viewer->getCamera();
    //osg::Vec3d vScreen(x,y, 0);
    osg::Matrix mVPW = camera->getViewMatrix() * camera->getProjectionMatrix() * camera->getViewport()->computeWindowMatrix();
    osg::Matrix invertVPW;
    invertVPW.invert(mVPW);
    vec3 = screenPoint * invertVPW;
    return vec3;
}

osg::Vec3d MyHandler::world2Camera(osg::Vec3 worldPoint, osgViewer::Viewer* viewer)//世界转相机
{
    osg::Vec3d vec3;
    osg::ref_ptr<osg::Camera> camera = viewer->getCamera();
    osg::Matrix mV = camera->getViewMatrix();
    vec3 = worldPoint * mV;
    return vec3;
}

osg::Vec3d MyHandler::camera2World(osg::Vec3 cameraPoint, osgViewer::Viewer* viewer)//相机转世界
{
    osg::Vec3d vec3;
    osg::ref_ptr<osg::Camera> camera = viewer->getCamera();
    //osg::Vec3d vScreen(x,y, 0);

    osg::Matrix mV = camera->getViewMatrix();
    osg::Matrix invertmV;
    invertmV.invert(mV);
    vec3 = cameraPoint * invertmV ;
    return vec3;
}
