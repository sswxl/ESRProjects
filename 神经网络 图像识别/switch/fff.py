import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# torch.load('filename.pth').to(device)

model = torch.load('343018.pkl', map_location=device)
model.eval()
batch_size = 1  #批处理大小
input_shape = (3,224,224)   #输入数据

input_data_shape = torch.randn(batch_size, *input_shape, device=device)

torch.onnx.export(model, input_data_shape, "tryyy.onnx", verbose=True)
