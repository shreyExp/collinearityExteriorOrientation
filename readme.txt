How to use the program

Requirements:
This program can run on any Windows and Linux operating system  with Python3 installed.
Python3 must have the numpy module installed.


Steps to run:
1. Prepare the points.json file:
This file will contain the GCP information. At least 3 are required. More points will qualify for
least aquared solution. A sample of the file is attached in the mail.

2. Prepare the Camera Parameters file:
This file (cameraParameters.json) will contain the focal legth of the camera and the pricipal point location on the image.
 A sample of this file is attached in the mail.

3. On a linux terminal type the following command:
python3 main.py > report.txt
This will run the main.py program. This program will automatically read the points.json and cameraParameters.json files.
After calculations. It will give the output in the report.txt file


