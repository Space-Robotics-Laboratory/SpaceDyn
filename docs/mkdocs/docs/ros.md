## 4. How to Execute Simulation
ClimbLabはMATLABで記述されており，実行にはMATALBを使用する必要がある．動作確認がおこなわれているMATLABのバージョンは，本論文執筆時(2022年2月1日)において，MATLAB R2019bである．よって，実行にはR2019b以降のMATLABがインストールされている環境が必要である．また，MATLAB Toolboxについては，基本的に必要としないが，ごく一部の機能について以下のToolbox：

* Optimization Toolbox
* Robotics Toolbox
* Signal Processing Toolbox
* Reinforcement Learning Toolbox
* Deep Learning Toolbox

が使用されるため，必要な場合は適宜インストールする必要がある．

### 4.1 How to execute `main_sim.m`
`main_sim.m`を実行し，歩行シミュレーションをおこなう手順について述べる．なお，先に述べたMATLABがインストールされている環境を前提とする．

まず基本設定である"default"で実行する場合は，

1. GitHubやBitBucketからClimbLabをダウンロードする．
2. MATLABを開く．
3. フォルダ/climblab/src下におかれている`main_sim.m`を開く．
4. `main_sim.m`内で使用するconfigが"default"に設定されていることを確認する．
5. MATLAB内の「現在のフォルダ」を`/climblab`にする．
6. `main\_sim.m`を実行する．
7. 実行後，「main_sim.mは現在のフォルダーやMATLABパス上で見つかりません．」というエラーが出た場合は，「パスに追加」を選択する．

という手順によって歩行シミュレーションが開始される．
<br>

またこの"default"から自分用に新たにシミュレーション設定を変更して試す場合には，

1. フォルダ`/climblab/config/USER`下におかれている`config\_USER\_param\_template.m`をコピーし，同じフォルダに，`config_USER_param.m`と名前を変えて作成する．
2. `config_USER_param.m`内に，defaultから変更したいシミュレーション設定を書き加える．
3. `main_sim.m`内で使用するconfigを"default"から"USER"に変更する．
4. `main_sim.m`を実行する．

という手順で設定を変更してシミュレーションをすることができる．変更できる設定は，フォルダ`config/default`に置かれている各classの初期設定用の関数の中身を見ることで確認できる．
<br>
また，自分が開発を進める中で，設定を半永続的に保存しておきたいconfigファイルがある場合は，presetとしてファイルを作成する．一方で，自分が開発したわけではない既存のpresetについては基本的に編集してはいけない．これらのpresetには論文などで対外的に示した結果を保存しておくべきものもあるので，編集する際は開発者やチームと話し合う必要がある．presetの作成手順を以下に示す．

1. `config/preset`下に，`config_xxx_param.m`という自分のconfigファイルを作る．**xxx**は任意だが，`config_`で始まり，`_param`で終わるように命名する．
2. `config_xxx_param.m`の内部に自分のパラメータを設定する．
3. main_sim内で
`config = ‘xxx’;
と設定する．
4. config_simulation内で一番下に，
`if strcmp(config, ‘xxx’)`
`[... ] = config_xxx_param(....);`
`end`
と新しく記述して，作成したconfigファイルをloadできるようにする．
5. `main_sim.m`内で使用するconfigを"xxx"に変更する．
6. `main_sim`を`/climblab`下で実行する．
