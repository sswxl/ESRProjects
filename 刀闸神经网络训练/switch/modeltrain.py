import torch as t
import torchvision as tv
import os
import matplotlib.pyplot as plt
import numpy as np
# 导入pytorch,numpy,matplotlib
from tqdm import tqdm


class DefaultConfigs(object):
    # 输入参数
    data_dir = "./imgs/"
    data_list = ["train", "test"]

    lr = 0.001
    epochs = 30
    num_classes = 3
    image_size = 224
    batch_size = 18
    channels = 3
    gpu = "0"

    use_gpu = t.cuda.is_available()


# 启用GPU

config = DefaultConfigs()
# 调用参数

normalize = tv.transforms.Normalize(mean=[0.485, 0.456, 0.406],
                                    std=[0.229, 0.224, 0.225]
                                    )
# 对图像进行规范化处理

transform = {
    config.data_list[0]: tv.transforms.Compose(
        [tv.transforms.Resize([224, 224]), tv.transforms.CenterCrop([224, 224]),
         tv.transforms.ToTensor(), normalize]  # 训练集的图像用CentrCrop进行剪裁
    ),
    config.data_list[1]: tv.transforms.Compose(
        [tv.transforms.Resize([224, 224]),
         tv.transforms.ToTensor(), normalize]  # 测试集图像不进行剪裁
    )
}
# 设置图像传入后的格式，转换为tensor的形式

datasets = {
    x: tv.datasets.ImageFolder(root=os.path.join(config.data_dir, x), transform=transform[x])
    for x in config.data_list
}  # 数据集

dataloader = {
    x: t.utils.data.DataLoader(dataset=datasets[x],
                               batch_size=config.batch_size,
                               shuffle=True
                               )
    for x in config.data_list
}  # 数据集加载器


def get_model(num_classes):
    # 建立神经网络
    model = tv.models.resnet34(pretrained=True)
    for parma in model.parameters():
        parma.requires_grad = False
    model.fc = t.nn.Sequential(
        t.nn.Dropout(p=0.3),
        t.nn.Linear(512, num_classes)
    )
    return (model)


def draw_lines(x, y1, y2, y3, y4):
    # 将训练所得的结果残差和准确率通过matplot画图输出
    y1 = np.array(y1)
    y2 = np.array(y2)
    y3 = np.array(y3)
    y4 = np.array(y4)

    # plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(x, y1, color='r', linestyle='-', marker="s", linewidth=1, label='train_loss')
    plt.plot(x, y3, color='b', linestyle='-', marker="s", linewidth=1, label='test_loss')
    plt.ylabel("loss")
    plt.title('Loss of model')
    plt.legend(['train_loss', 'test_loss'], loc='upper right')

    # plt.figure(2)
    plt.subplot(2, 1, 2)
    plt.plot(x, y2, color='r', linestyle='-', marker="x", linewidth=1, label='train_acc')
    plt.plot(x, y4, color='b', linestyle='-', marker="x", linewidth=1, label='test_acc')
    plt.xlabel("epcho")
    plt.ylabel("acc")
    plt.title('Accuracy of model')
    plt.legend(['train_acc', 'test_acc'], loc='lower right')

    plt.show()


def train(epochs):
    # 训练模型
    model = get_model(config.num_classes)  # 调用模型
    print(model)
    loss_f = t.nn.CrossEntropyLoss()
    if (config.use_gpu):  # 如果gpu能用调用gpu
        model = model.cuda()
        loss_f = loss_f.cuda()

    opt = t.optim.Adam(model.fc.parameters(), lr=config.lr)  # 调用加速器

    acctrain = 0
    tr_loss = []
    tr_acc = []
    te_loss = []
    te_acc = []
    for epoch in tqdm(range(epochs)):
        train_loss = []
        train_acc = []
        test_loss = []
        test_acc = []
        model.train(True)  # 模式设置为训练模式，对模型进行训练
        print("Epoch {}/{}".format(epoch + 1, epochs))
        for batch, datas in enumerate(iter(dataloader["train"])):
            x, y = datas
            if (config.use_gpu):
                x, y = x.cuda(), y.cuda()
            y_ = model(x)
            # print(x.shape,y.shape,y_.shape)
            _, pre_y_ = t.max(y_, 1)
            pre_y = y
            # print(y_.shape)
            loss = loss_f(y_, pre_y)
            # print(y_.shape)
            acc = t.sum(pre_y_ == pre_y)

            loss.backward()
            opt.step()
            opt.zero_grad()
            if (config.use_gpu):
                loss = loss.cpu()
                acc = acc.cpu()
            train_loss.append(loss.data)
            train_acc.append(acc)
            # if((batch+1)%5 ==0):

        print("Batch {}, Train loss:{:.4f}, Train acc:{:.4f}"
              .format(batch + 1, np.mean(train_loss) / config.batch_size, np.mean(train_acc) / config.batch_size))
        tr_acc.append(np.mean(train_acc) / config.batch_size)
        tr_loss.append(np.mean(train_loss) / config.batch_size)

        if acctrain < np.mean(train_acc) / config.batch_size:
            acctrain = np.mean(train_acc) / config.batch_size
            print(acctrain)
            dowloadmodle = True
        else:
            dowloadmodle = False

        model.train(False)  # 关闭训练模式，调用测试集进行测试

        for batch, datas in enumerate(iter(dataloader["test"])):
            x, y = datas
            if (config.use_gpu):
                x, y = x.cuda(), y.cuda()
            y_ = model(x)
            # print(x.shape,y.shape,y_.shape)
            _, pre_y_ = t.max(y_, 1)
            pre_y = y
            # print(y_.shape)
            loss = loss_f(y_, pre_y)
            acc = t.sum(pre_y_ == pre_y)

            if (config.use_gpu):
                loss = loss.cpu()
                acc = acc.cpu()

            test_loss.append(loss.data)
            test_acc.append(acc)

        print("Batch {}, Test loss:{:.4f}, Test acc:{:.4f}"
              .format(batch + 1, np.mean(test_loss) / config.batch_size, np.mean(test_acc) / config.batch_size))
        te_acc.append(np.mean(test_acc) / config.batch_size)
        te_loss.append(np.mean(test_loss) / config.batch_size)

        if dowloadmodle:
            t.save(model, str(epoch + 1) + "resnet34.pkl")

    x = np.arange(0, epoch + 1, 1)
    draw_lines(x, tr_loss, tr_acc, te_loss, te_acc)


# 训练开始
if __name__ == "__main__":
    train(config.epochs)
