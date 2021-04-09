import torch
import torchvision
from torchvision import transforms
from PIL import Image
from unis import *
import time

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

win_name='detect'

# download model
model=torch.load('503018.pkl')
# model=torchvision.models.mobilenet_v2(pretrained=True)

img_size=224

normalize=transforms.Normalize(mean=[0.485,0.456,0.406],
                                  std=[0.229,0.224,0.225])

transform=transforms.Compose(
    [transforms.Resize([img_size,img_size]),
     transforms.CenterCrop([img_size,img_size]),
     transforms.ToTensor(),
     normalize]
)

model=model.eval()
model=model.to(device)
cam = cv2.VideoCapture(0)

classes=['close','full-open','half-open']

with torch.no_grad():
    while True:
        ret, original = cam.read()
        #cv2.imshow(win_name, original)
        recent_time=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
        mat = cv2.cvtColor(original, cv2.COLOR_BGR2RGB)
        mat_test = Image.fromarray(mat)
        frame = transform(mat_test).unsqueeze(0)
        frame=frame.to(device)
        prediction=model(frame)
        try:
            _, indices = torch.sort(prediction, descending=True)
            percentage = torch.nn.functional.softmax(prediction, dim=1)[0] * 100
            item = indices[0]

            _, predicted = torch.max(prediction.data, 1)
            maxPredicted = classes[predicted.item()]
            maxAccuracy = percentage[item[0]].item()

            ShowPicResult(original, win_name, maxAccuracy, maxPredicted, recent_time, False)
        except:
            print('something went wrong!')
        if cv2.waitKey(20) & 0xFF == ord('q'):
            break

    cam.release()
    cv2.destroyAllWindows()
