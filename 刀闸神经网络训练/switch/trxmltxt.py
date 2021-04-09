import os
import os.path
import glob
# xml为自带库
import xml.etree.ElementTree as ET


# 训练集信息
class_names = ['switch close', 'switch full open','switch half open']	# 训练集标签, 比如 dog, cat 

# 获取文件路径
# xml, txt本地路径
xml_path = './test labels/'
txt_path = './txts/'

def xml_to_txt_onefiles(xml_path, txt_path):
    # 创建txt文件
    name_txt = 'test.txt'
    txt_file = os.path.join(txt_path, name_txt)
    f = open(txt_file, 'w')

    # 寻找xml文件, 为空抛出异常 FileExistsError
    try:
        annotations = os.listdir(xml_path)
        annotations = glob.glob(str(xml_path) + "\\" + str(annotations) + '*.xml')

        if not annotations:
            raise FileExistsError
    except FileExistsError:
        print("No xml files, check the path!")

    # 遍历xml文件
    for _, file in enumerate(annotations):
        xml_file = open(file)
        tree = ET.parse(xml_file)	# 读取xml文件
        root = tree.getroot()		# 读取root节点

        # 寻找所需信息
        filename = root.find('filename').text
        filepath = root.find('path').text
        f.write(filepath)
        for obj in root.iter('object'):
            name = obj.find('name').text
            class_num = class_names.index(name)

            bndbox = obj.find('bndbox')
            x1 = bndbox.find('xmin').text
            x2 = bndbox.find('xmax').text
            y1 = bndbox.find('ymin').text
            y2 = bndbox.find('ymax').text

            # 写入和关闭txt
            f.write(' '+x1+','+y1+','+x2+','+y2+','+str(class_num))
        f.write('\n')

if __name__ == '__main__':
    xml_to_txt_onefiles(xml_path, txt_path)
    print('Finish!!')
