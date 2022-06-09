# This Python 3 environment comes with many helpful analytics libraries installed
# It is defined by the kaggle/python Docker image: https://github.com/kaggle/docker-python
# For example, here's several helpful packages to load

import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import plotly.graph_objects as go
import plotly.offline as py

# Input data files are available in the read-only "../input/" directory
# For example, running this (by clicking run or pressing Shift+Enter) will list all files under the input directory

import os
for dirname, _, filenames in os.walk('/kaggle/input'):
    for filename in filenames:
        print(os.path.join(dirname, filename))

# You can write up to 5GB to the current directory (/kaggle/working/) that gets preserved as output when you create a version using "Save & Run All" 
# You can also write temporary files to /kaggle/temp/, but they won't be saved outside of the current session
df = pd.read_excel('/kaggle/input/komp-plasma-metabolomics-dataset/Steroids_Bile_Acid_details.xlsx')
df.head()
df.isnull().sum()
plt.figure(figsize = (10,6))
sns.lineplot(data=df, x="PubChem CID", y="CXP Collision Cell Exil Potential", color='r', marker='*')
plt.title("CXP Collision Cell Exil Potential CID",fontweight='bold',size=20)
plt.xlabel('PubChem CID',size=15)
plt.ylabel('CXP Collision Cell Exil Potential',size=15)
plt.xticks(rotation=45)
plt.show()

plt.figure(figsize = (10,6))
sns.lineplot(data=df, x="PubChem CID", y="EP Entrance Potential", color='b', marker='*')
plt.title("EP Entrance Potential CID",fontweight='bold',size=20)
plt.xlabel('PubChem CID',size=15)
plt.ylabel('EP Entrance Potential',size=15)
plt.xticks(rotation=45)
plt.show()

#Let's visualise steroids
steroids = df.groupby('PubChem CID').sum()[['DP', 'CE', 'EP Entrance Potential']]
#evolution['Expiration Rate'] = (evolution['Expired'] / evolution['Cumulative']) * 100
#evolution['Discharging Rate'] = (evolution['Discharged'] / evolution['Cumulative']) * 100
steroids.head()

#DP: Declustering Potential (V); EP: Entrance Potential (V); CE: Collision Energy

plt.figure(figsize=(20,7))
plt.plot(steroids['CE'], label='Collision Energy')
plt.plot(steroids['DP'], label='Declustering Potential')
plt.plot(steroids['EP Entrance Potential'], label='EP Entrance Potential')
plt.legend()
#plt.grid()
plt.title('Collision Energy, Declustering & Entrance Potentials')
plt.xticks(steroids.index,rotation=45)
plt.xlabel('PubChem CID')
plt.ylabel('Count')
plt.show()

#What about disaggregated
plt.figure(figsize=(20,7))
plt.plot(steroids['DP'], label='EP Entrance Potential')
plt.legend()
#plt.grid()
plt.title('EP Entrance Potential')
plt.xticks(steroids.index,rotation=45)
plt.ylabel('Count')
plt.show()

#This is another way of visualizing the sex-disaggregated data
diff_steroids = steroids.diff().iloc[1:]
plt.figure(figsize=(20,7))
plt.plot(diff_steroids['DP'], label='Declustering Potential')
plt.legend()
plt.grid()
plt.title('Declustering Potential')
plt.xticks(steroids.index,rotation=45)
plt.ylabel('Count')
plt.show()

plt.figure(figsize=(20, 10))
plt.subplot(431)
sns.countplot(df['ESI mode'])
plt.title('Electrospray Ionization')
plt.xlabel('')
plt.subplot(432)
sns.countplot(df['MRM'])
plt.title('Multiple Reactions Monitoring')
plt.xlabel('MRM')
plt.xticks(rotation=45)
plt.subplot(433)
sns.countplot(df['LOD (nM)'])
plt.title('Limit of Detection (nM)')
plt.xlabel('LOD (nM)')
plt.xticks(rotation=45)
plt.subplot(434)
sns.countplot(df['LOQ (nM)'])
plt.title('Limit of Quantification (nM)')
plt.xlabel('Limit of Quantification (nM)')

plt.figure(figsize=(20, 12))
plt.subplot(341)
sns.boxplot(x=df['DP'])
plt.title('Declustering Potencial')
plt.xlabel('DP')
plt.subplot(342)
sns.boxplot(x=df['CE'])
plt.title('Collision Energy')
plt.xlabel('CE')
plt.subplot(343)
sns.boxplot(x=df['EP Entrance Potential'])
plt.title('EP Entrance Potential')
plt.xlabel('EP Entrance Potential')
plt.subplot(344)
sns.boxplot(x=df['CXP Collision Cell Exil Potential'])
plt.title('Cell Exit Potential')
plt.xlabel('CXP Collision Cell Exil Potential')

plt.figure(figsize=(20,8))
plt.subplot(211)
sns.distplot(df['acc.mass (neutral)'])
plt.title('Accurate Mass Neutral')
plt.xlabel('')
plt.subplot(212)
sns.distplot(df['acc mass [M+H]+ resp. [M-H]-'])
plt.title('Accurate Mass M+H ions')
plt.xlabel('acc mass [M+H]+ resp. [M-H]-')

#word cloud
from wordcloud import WordCloud, ImageColorGenerator
text = " ".join(str(each) for each in df.name)
# Create and generate a word cloud image:
wordcloud = WordCloud(max_words=200,colormap='Set3', background_color="black").generate(text)
plt.figure(figsize=(10,6))
plt.figure(figsize=(15,10))
# Display the generated image:
plt.imshow(wordcloud, interpolation='Bilinear')
plt.axis("off")
plt.figure(1,figsize=(12, 12))
plt.show()