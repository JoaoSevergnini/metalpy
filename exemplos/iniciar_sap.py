import sys
import comtypes.client

def IniciarSAP():

    AttachToInstance = False

    SpecifyPath = True
    
    ProgramPath = 'C:\Program Files\Computers and Structures\SAP2000 24\SAP2000.exe'


    if AttachToInstance:

        # attach to a running instance of SAP2000

        try:

            # get the active SapObject

            mySapObject = comtypes.client.GetActiveObject("CSI.SAP2000.API.SapObject")

        except (OSError, comtypes.COMError):

            print("No running instance of the program found or failed to attach.")

            sys.exit(-1)

    else:

        # create API helper object

        helper = comtypes.client.CreateObject('SAP2000v1.Helper')

        helper = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)

        if SpecifyPath:

            try:

                # 'create an instance of the SAPObject from the specified path

                mySapObject = helper.CreateObject(ProgramPath)

            except (OSError, comtypes.COMError):

                print("Cannot start a new instance of the program from " + ProgramPath)

                sys.exit(-1)

        else:

            try:

                # create an instance of the SAPObject from the latest installed SAP2000

                mySapObject = helper.CreateObjectProgID("CSI.SAP2000.API.SapObject")

            except (OSError, comtypes.COMError):

                print("Cannot start a new instance of the program.")

                sys.exit(-1)

        # start SAP2000 application

        mySapObject.ApplicationStart()

    # create SapModel object

    SapModel = mySapObject.SapModel


    return SapModel