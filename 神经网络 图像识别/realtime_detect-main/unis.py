import cv2
import tkinter as tk
import numpy as np

#def GetScreenCenter():
    #root = tk.Tk()
    #return root.winfo_screenwidth()//2,root.winfo_screenheight()//2

#def AdaptSize(img):
    # 视频、图片过大直接1/2
    #center_x, center_y = GetScreenCenter()
    #img_h, img_w, _ = img.shape
    #if img_h > center_y * 2 or img_w > center_x * 2:
    #    img = cv2.resize(img, (img_w // 2, img_h // 2))
    #return img

def CentralShow(win_name,img,stop):
    #center_x, center_y = GetScreenCenter()
    #img=AdaptSize(img)
    #img_h,img_w,_=img.shape
    #t_x, t_y = (center_x - img_w // 2), (center_y - img_h // 2)
    #cv2.namedWindow('picture', 0)
    cv2.imshow(win_name, img)
    #cv2.moveWindow(win_name, t_x, t_y)
    if stop:
        cv2.waitKey(0)

def ShowPicResult (cv_mat, win_name, maxpercentage, labels, r_t, stop = False):
    cv_mat = cv2.putText(cv_mat, labels, (0, 30), cv2.FONT_HERSHEY_TRIPLEX, 1, (0, 255, 0), 1)
    cv_mat = cv2.putText(cv_mat, 'Accuracy'+str(maxpercentage), (0, 60), cv2.FONT_HERSHEY_TRIPLEX, 0.7, (0, 255, 0), 1)
    cv_mat = cv2.putText(cv_mat, r_t, (450, 15), cv2.FONT_HERSHEY_TRIPLEX, 0.5, (0, 255, 0), 1)
    CentralShow(win_name, cv_mat, stop)
