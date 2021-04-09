import torch
import torchvision
from torchvision import transforms
from PIL import Image
from unis import *

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")#调用gpu

win_name='detect'
classes=['close','full-open','half-open']

model= torch.load('343018.pkl')#载入模型

img_size=224

normalize=transforms.Normalize(mean=[0.485,0.456,0.406],
                                  std=[0.229,0.224,0.225])

transform=transforms.Compose(
    [transforms.Resize([img_size,img_size]),
     transforms.CenterCrop([img_size,img_size]),
     transforms.ToTensor(),
     normalize]
)

# change the path to the one you need
mat=cv2.imread('IMG_20201221_100219.jpg')
mat=cv2.cvtColor(mat,cv2.COLOR_BGR2RGB)
mat_test=Image.fromarray(mat)

test=transform(mat_test).unsqueeze(0)

model=model.to(device)
test=test.to(device)

with torch.no_grad():
    predict=model(test)

_,indices=torch.sort(predict,descending=True)
percentage=torch.nn.functional.softmax(predict,dim=1)[0]*100
item=indices[0]

_,predicted=torch.max(predict.data,1)
maxPredicted = classes[predicted.item()]
maxAccuracy = percentage[item[0]].item()
print('start')
ShowPicResult(mat,win_name,maxAccuracy,maxPredicted,True)


