import torch as t
import torchvision as tv
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
#调用模块

device=t.device('cuda:0' if t.cuda.is_available() else "cpu")

img_size=224

normalize=tv.transforms.Normalize(mean=[0.485,0.456,0.406],
                                  std=[0.229,0.224,0.225])

transform=tv.transforms.Compose(
    [tv.transforms.Resize([img_size,img_size]),tv.transforms.CenterCrop([img_size,img_size]),
     tv.transforms.ToTensor(),
     normalize]
)

path=data_dir="./imgs/test/"
testset=tv.datasets.ImageFolder(path,transform=transform)
testloader=t.utils.data.DataLoader(testset,batch_size=40,shuffle=True,num_workers=6)

#classes=['close','half-open','full-open']
classes=['0','1','2']
type_cl=['close','full-open','half-open']

#def imshow(img):
    #plt.imshow(img)
    #plt.show()

def predic(model,img_path,imgType,isShowSoftmax=False,isShowImg=False):
    #预测单个图片十分正确
    t.no_grad()
    image_PIL=Image.open(img_path)

    image_tensor=transform(image_PIL)
    image_tensor=t.unsqueeze(image_tensor,dim=0).float()

    #image_tensor=image_tensor.cuda()
    out=model(image_tensor)
    print(out)

    _,indices=t.sort(out,descending=True)
    percentage = t.nn.functional.softmax(out, dim=1)[0] * 100
    print(percentage)

    item=indices[0]
    if isShowSoftmax:
        for idx in item:
            ss=percentage[idx]
            value=ss.item();
            name=classes[idx]
            print('name:',name,'predict percentage:',value)

    _, predicted = t.max(out.data, 1)
    print(predicted)
    maxPredicted = classes[predicted.item()]
    print(maxPredicted)
    maxAccuracy = percentage[item[0]].item()
    print(maxAccuracy)
    if imgType==maxPredicted:
        if maxPredicted=='0':
            predict='close'
        elif maxPredicted == '1':
            predict = 'fullopen'
        elif maxPredicted == '2':
            predict = 'halfopen'
        print('right predict is:',predict,'predict percentage',maxAccuracy)
    else:
        if maxPredicted=='0':
            predict='close'
        elif maxPredicted == '1':
            predict = 'fullopen'
        elif maxPredicted == '2':
            predict = 'halfopen'
        if imgType =='0':
            real='close'
        elif imgType == '1':
            real = 'fullopen'
        elif imgType == '2':
            real = 'halfopen'
        print('wrong predict is:',predict,'right is',real,'predict percentage',maxAccuracy)
    if isShowImg:
        plt.imshow(image_PIL)
        plt.show()

#def get_model(num_classes):
    #神经网络模型
#    model=tv.models.resnet50(pretrained=True)

#    model.fc=t.nn.Sequential(
#        t.nn.Dropout(p=0.3),
#        t.nn.Linear(2048,num_classes)
#    )
#    return model

def testAll(model):
    #从测试集中随机选择图片进行预测
    dataiter=iter(testloader)

    images,labels=dataiter.next()
    #if (t.cuda.is_available()):
        #images, labels = images.cuda(), labels.cuda()
    print(labels)
#    print('real:',' '.join('%5s'%type_cl[labels[j]]for j in range(40)))
    outputs=model(images)
    print(outputs)
    _,predicted=t.max(outputs.data,1)

#    print('predict:',' '.join('%5s'%type_cl[predicted[j]]for j in range(40)))

    right_num=0    
    for i in range(40):
        print(i+1)
        if type_cl[labels[i]]==type_cl[predicted[i]]:
            right_num+=1
            print('predict is right, the switch is:',type_cl[predicted[i]])
        else:
            print('predict is wrong, the predict is:',type_cl[predicted[i]],'but the real is:',type_cl[labels[i]])
    print('we have test 40 pic, and there are',str(right_num),'right print in this')

    imshow2(tv.utils.make_grid(images,nrow=5))

def imshow2(img):
    #展示图片
    #if (t.cuda.is_available()):
        #img=img.cpu()
    img=img/2+0.5
    npimg=img.numpy()
    plt.imshow(np.transpose(npimg,(1,2,0)))
    plt.show()


if __name__ == "__main__":
    #载入已经训练好的pkl文件
    model=t.load('./503018.pkl',map_location='cpu')

    #model=model.to(device)
    model.eval()

    #对单个图片进行预测
    predic(model,'imgs/test/close/my1.jpg', classes[0], False, True)
    predic(model,'imgs/test/close/IMG_20201221_102830.jpg', classes[0], False, True)
    predic(model,'imgs/test/close/IMG_20201221_102823.jpg', classes[0], False, True)
    predic(model,'imgs/test/fullopen/IMG_20201221_100432.jpg', classes[1], False, True)
    predic(model,'imgs/test/fullopen/IMG_20201221_103601.jpg', classes[1], False, True)
    predic(model,'imgs/test/halfopen/IMG_20201221_103737.jpg', classes[2], False, True)
    predic(model,'imgs/test/halfopen/IMG_20201221_153218.jpg', classes[2], False, True)
    predic(model,'imgs/train/fullopen/IMG_20201221_145940_1.jpg', classes[1], False, True)
    predic(model,'imgs/train/close/IMG_20201221_143037.jpg', classes[0], False, True)

    testAll(model)
