//==============================================================================
/*
    Simple haptic interface handler.
    Tristan Croll
    University of Cambridge
    2016
*/
//==============================================================================


//==============================================================================
/*
    Software License Agreement (BSD License)
    Copyright (c) 2003-2016, CHAI3D.
    (www.chai3d.org)

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

    * Neither the name of CHAI3D nor the names of its contributors may
    be used to endorse or promote products derived from this software
    without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
    ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

    \author    <http://www.chai3d.org>
    \author    Francois Conti
    \version   3.1.1 $Rev: 1869 $
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"

//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------



extern "C"
{
    //------------------------------------------------------------------------------
    // DECLARED FUNCTIONS
    //------------------------------------------------------------------------------

    // main haptics simulation loop
    void updateHaptics(void);

    // get the number of connected devices
    int getNumDevices(void);

    // set the spring constant for device i
    void setSpringConstant (int i, double k);

    // switch damping on/off for device i
    void setDamping(int i, bool d);

    // set the target position for device i
    void setTargetPosition(int i, double x, double y, double z);

    // Attach device i to an object
    void startTugging(int i);

    // Detach device i
    void stopTugging(int i);

    // Turn on force feedback for device i
    void turnOnFeedback(int i);

    // Turn off force feedback for device i
    void turnOffFeedback(int i);

    // return the (x,y,z) position of device i
    double* getPosition(int i);

    // get the states of the buttons on device i
    bool* getButtonStates (int i);



    //------------------------------------------------------------------------------
    // GENERAL SETTINGS
    //------------------------------------------------------------------------------

    // maximum number of devices supported by this application
    static const int MAX_DEVICES = 16;

    //------------------------------------------------------------------------------
    // DECLARED VARIABLES
    //------------------------------------------------------------------------------

    // a haptic device handler
    cHapticDeviceHandler* handler;

    // a pointer to the current haptic device
    cGenericHapticDevicePtr hapticDevice[MAX_DEVICES];

    // number of haptic devices detected
    int numHapticDevices = 0;

    // variable to store the position [m] of each haptic device
    cVector3d hapticDevicePosition[MAX_DEVICES];

    // variable to store the target position of each haptic device
    cVector3d targetPosition[MAX_DEVICES];

    // spring constants connecting each device to an object
    double springConstant[MAX_DEVICES] = {0};

    // struct holding information about each device
    cHapticDeviceInfo info[MAX_DEVICES];

    // 2D vector holding the states of four buttons for each device
    //vector<vector<bool>> hapticDeviceButtons(MAX_DEVICES, vector <bool>(4, false));

    // flag for using damping (ON/OFF)
    bool useDamping[MAX_DEVICES] = {false};

    // flag for using force feedback (ON/OFF)
    bool useFeedback[MAX_DEVICES] = {true};

    // is this device currently pulling something?
    bool deviceInUse[MAX_DEVICES] = {false};

    // flag to indicate if the haptic simulation currently running
    bool simulationRunning = false;

    // flag to indicate if the haptic simulation has terminated
    bool simulationFinished = true;

    // frequency counter to measure the simulation haptic rate
    cFrequencyCounter frequencyCounter;

    // main thread running a loop driving the haptic interface(s)
    cThread* hapticsThread;


    cHapticDeviceHandler* HapticHandler(void)
    {
        if (!handler)
        {
            //--------------------------------------------------------------------------
            // HAPTIC DEVICES
            //--------------------------------------------------------------------------

            // create a haptic device handler
            handler = new cHapticDeviceHandler();

            // get number of haptic devices
            numHapticDevices = handler->getNumDevices();

            // setup each haptic device
            for (int i=0; i<numHapticDevices; i++)
            {
                // get a handle to the first haptic device
                handler->getDevice(hapticDevice[i], i);

                // open a connection to haptic device
                hapticDevice[i]->open();

                // calibrate device (if necessary)
                hapticDevice[i]->calibrate();

                // retrieve information about the current haptic device
                info[i] = hapticDevice[i]->getSpecifications();



                // display a reference frame if haptic device supports orientations
                if (info[i].m_sensedRotation == true)
                {
                    // what can we do with a haptic device that allows rotations?
                }

                // if the device has a gripper, enable the gripper to simulate a user switch
                hapticDevice[i]->setEnableGripperUserSwitch(true);

                // initialise the target position
                targetPosition[i].set(0.0, 0.0, 0.0);

            }


            //--------------------------------------------------------------------------
            // START SIMULATION
            //--------------------------------------------------------------------------
            // create a thread which starts the main haptics rendering loop
            hapticsThread = new cThread();
            hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);
        }
        return handler;

    }



    //------------------------------------------------------------------------------

    void stopHaptics(void)
    {
        // stop the simulation
        simulationRunning = false;

        // wait for haptics loop to terminate
        while (!simulationFinished) { cSleepMs(100); }

        // close haptic devices
        for (int i=0; i<numHapticDevices; i++)
        {
            hapticDevice[i]->close();
        }

        delete hapticsThread;
        delete handler;
    }


    //------------------------------------------------------------------------------

    // get the number of connected devices
    int getNumDevices(void)
    {
        return numHapticDevices;
    }

    // set the spring constant for device i
    void setSpringConstant (int i, double k)
    {
        springConstant[i] = k;
    }

    // switch damping on/off for device i
    void setDamping(int i, bool d)
    {
        useDamping[i] = d;
    }

    // set the target position for device i
    void setTargetPosition(int i, double x, double y, double z)
    {
        targetPosition[i].set(x, y, z);
    }

    // Attach device i to an object
    void startTugging(int i)
    {
        deviceInUse[i] = true;
    }

    // Detach device i
    void stopTugging(int i)
    {
        deviceInUse[i] = false;
    }

    // Turn on force feedback for device i
    void turnOnFeedback(int i)
    {
        useFeedback[i] = true;
    }

    // Turn off force feedback for device i
    void turnOffFeedback(int i)
    {
        useFeedback[i] = false;
    }

    // return the (x,y,z) position of device i
    double* getPosition(int i)
    {
        double* pos = new double[3];
        pos[0] = hapticDevicePosition[i].x();
        pos[1] = hapticDevicePosition[i].y();
        pos[2] = hapticDevicePosition[i].z();

        return pos;
    }

    // Return the states of up to four buttons on the haptic device
    bool* getButtonStates(int i)
    {
        bool * buttons = new bool[4];

        for (int j=0; j < 4; j++) {
            hapticDevice[i]->getUserSwitch(j, buttons[j]);
        }

        return buttons;

    }

    // Main loop run by cThread object to update haptic device position(s) and
    // manage force feedback
    void updateHaptics(void)
    {
        // initialize frequency counter
        frequencyCounter.reset();

        // simulation in now running
        simulationRunning  = true;
        simulationFinished = false;

        // main haptic simulation loop
        while(simulationRunning)
        {
            for (int i=0; i<numHapticDevices; i++)
            {
                /////////////////////////////////////////////////////////////////////
                // READ HAPTIC DEVICE
                /////////////////////////////////////////////////////////////////////

                // read position
                cVector3d position;
                hapticDevice[i]->getPosition(position);

                // read orientation
                cMatrix3d rotation;
                hapticDevice[i]->getRotation(rotation);

                // read linear velocity
                cVector3d linearVelocity;
                hapticDevice[i]->getLinearVelocity(linearVelocity);

                // read angular velocity
                cVector3d angularVelocity;
                hapticDevice[i]->getAngularVelocity(angularVelocity);

                // update global variable for graphic display update
                hapticDevicePosition[i] = position;

                // If this device isn't currently dragging something, there's no need to
                // compute forces for it
                if (!deviceInUse[i]) { continue; }
                /////////////////////////////////////////////////////////////////////
                // COMPUTE AND APPLY FORCES
                /////////////////////////////////////////////////////////////////////


                // desired orientation
                cMatrix3d desiredRotation;
                desiredRotation.identity();

                // variables for forces
                cVector3d force (0,0,0);
                cVector3d torque (0,0,0);

                bool canRotate = info[i].m_sensedRotation;


                // apply force field
                if (useFeedback[i])
                {
                    cVector3d desiredPosition = targetPosition[i];
                    // compute linear force
                    double Kp = springConstant[i]; // [N/m]
                    cVector3d forceField = Kp * (desiredPosition - position);
                    force.add(forceField);

                    // apply damping term
                    if (useDamping)
                    {
                        // compute linear damping force
                        double Kv = 1.0 * info[i].m_maxLinearDamping;
                        cVector3d forceDamping = -Kv * linearVelocity;
                        force.add(forceDamping);
                    }

                    if (canRotate)
                    {
                        // compute angular torque
                        double Kr = 0.05; // [N/m.rad]
                        cVector3d axis;
                        double angle;
                        cMatrix3d deltaRotation = cTranspose(rotation) * desiredRotation;
                        deltaRotation.toAxisAngle(axis, angle);
                        torque = rotation * ((Kr * angle) * axis);

                        if (useDamping)
                        {
                            // compute angular damping force
                            double Kvr = 1.0 * info[i].m_maxAngularDamping;
                            cVector3d torqueDamping = -Kvr * angularVelocity;
                            torque.add(torqueDamping);
                        }

                    }
                }


                // send computed force, torque, and gripper force to haptic device
                hapticDevice[i]->setForceAndTorque(force, torque);
            }

            // update frequency counter
            frequencyCounter.signal(1);
        }

        // exit haptics thread
        simulationFinished = true;
    }

}
