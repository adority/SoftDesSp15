""" Experiment with face detection and image filtering using OpenCV """

import cv2
import numpy as np

face_cascade = cv2.CascadeClassifier('/home/aubrey/Documents/SoftDes/SoftDesSp15/toolbox/image_processing/haarcascade_frontalface_alt.xml')
kernel = np.ones((10,10),'uint8')
cap = cv2.VideoCapture(0)

while(True):
	ret, frame = cap.read()

	faces = face_cascade.detectMultiScale(frame,scaleFactor=1.2, minSize=(20,20))
	for (x,y,w,h) in faces:
		frame[y:y+h, x:x+w, :] = cv2.dilate(frame[y:y+h, x:x+w, :], kernel)
		cv2.circle(frame, (int(x+(.3*w)),int(y+(.3*h))), int(.1*w), (255,255,255), -1, 8, 0)
		cv2.circle(frame, (int(x+(.7*w)),int(y+(.3*h))), int(.1*w), (255,255,255), -1, 8, 0)
		cv2.circle(frame, (int(x+(.3*w)),int(y+(.3*h))), int(.05*w), (0,0,0), -1, 8, 0)
		cv2.circle(frame, (int(x+(.7*w)),int(y+(.3*h))), int(.05*w), (0,0,0), -1, 8, 0)
		cv2.ellipse(frame, (int(x+(.5*w)),int(y+(.6*h))), (int(.4*w),int(.3*h)),0,0,180,(0,0,0), -1, 8, 0)
		cv2.rectangle(frame, (int(x+(.3*w)),int(y+(.6*h))), (int(x+(.4*w)),int(y+(.7*h))), (255,255,255), -1)
		#cv2.rectangle(frame,(x,y),(x+w,y+h),(0,0,255))

	cv2.imshow('frame', frame)

	if cv2.waitKey(1) & 0xFF == ord('q'):
		break

cap.release()
cv2.destroyAllWindows()