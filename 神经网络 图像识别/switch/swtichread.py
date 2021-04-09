import torch as t
import torchvision as tv
import matplotlib.pyplot as plt
from PIL import Image
from mynn import *
import cv2

model = t.load('343018.pkl').cuda()  #加载训练好的pytorch模型

normalize = tv.transforms.Normalize(mean = [0.485, 0.456, 0.406],
                                    std = [0.229, 0.224, 0.225]
                                    )

transforms = tv.transforms.Compose([
    tv.transforms.Grayscale(),
    tv.transforms.Resize([224,224]),
    tv.transforms.ToTensor(),
    normalize
])#定义图像变换以符合网络输入

switch=['close','full-open','half-open']

cap = cv2.VideoCapture(0)# 摄像头，0是笔记本自带摄像头

while(cap.isOpened()):
    ret,frame=cap.read()
    frame=frame[:,::-1,:]
    frame=frame.copy()
    gray=cv2.cvtColor(frame,cv2.COLOR_BGR2GRAY)
    img=frame
    img=Image.fromarray(img)
    img=transforms(img)
    pre=model(img).max(1)[1].item()
    frame=cv2.putText(frame,switch[pre],(100,100),cv2.FONT_HERSHEY_SIMPLEX,1,(55,225,155),2)
    cv2.imshow('switch',frame)
    if cv2.waikey(1)==ord('q'):
        break
cap.release()
cv2.destroyWindows()
