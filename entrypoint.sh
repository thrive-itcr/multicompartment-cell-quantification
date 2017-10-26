#!/bin/sh
# Copyright (c) General Electric Company, 2017.  All rights reserved.

/usr/bin/python ./rt106GenericAdaptorREST.py & sleep 3

/usr/bin/python ./rt106GenericAdaptorAMQP.py --broker rabbitmq --dicom http://datastore:5106
