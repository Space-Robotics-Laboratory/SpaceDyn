## 5. Folder Tree
本章では，ClimbLabのフォルダ構成と関数の置き場所を示す．なお，太字はフォルダ名を表す．本章を軽く一読することで，ClimbLab全体の構造を把握することできる．また，新たな関数やファイルを作成する場合には，このフォルダ構成を参考に適切な場所に配置すること．

** climblab/ **

- config config関数群一式を保存するフォルダ

	- **default** : configのdefault用関数群
  
	- **preset** : これまでの開発に使われたconfiguration関数群
  
	- **USER** : 変数上書き用USER config関数とそのテンプレート
  
	- `config_simulation.m` : シミュレーションのconfigを設定する大元の関数
  
	- `config_all_default_param.m` : defaultのなかのconfigファイル群を実行する関数

- **dat** : dataフォルダーの意．実行したシミュレーションのログや動画の保存先

- **docs** 

	- media : 対外用に使用する画像・gifの保存場所
  
	- mkdocs : 本マニュアルの保存場所
  
- **lib** : ClimbLabで使用している関数の中で，脚型ロボットの歩行のみならず他の分野への汎用性も高い関数群を保存する場所

	- **equ_gia** : GIA安定領域についての計算に関わる関数一式
  
	- **spacedyn_v2r0** : SpaceDynの関数一式
  
	- **MATLAB Reference** : 各関数の役割について，html形式で見やすくまとめられたもの
  
    - target_detection : 把持候補点検出に関わる関数一式
  
    - 回転行列を計算するだけの関数群 : `rot_2xy.m`, `rot_x.m`, `rot_y.m, rot_z.m`

- **src** : 脚型クライミングロボットの歩行シミュレーションに特化して開発されてきた関数群．下位の構成フォルダの数が多いので，src以下にのみ項目をわけて，次に改めて記す．

<br>

**climblab/src/**

- **controller** : ループ内で使用される歩容計画に関わる重要な関数の多くが置かれている

	- **control** : 関節トルクを制御するための関数を収容
  
	- **equilibrium** : 安定性解析に関わる関数を収容．
  
	- **gait** : 歩容計画に関わる一連の関数群を収容．実行される順番に従って下位フォルダに番号が振られている．各フォルダ内に，各歩容計画タイプに応じた関数ファイルが入っている．
  
		- **0_upd_graspable_points_in_reachable_area** : 可動範囲内の把持候補点の位置取得
  
		- **1_upd_swing_num** : 遊脚の選定をする関数
  
		- **2_upd_swing_next_pos** : 次の遊脚の接地点(把持点)を決定する関数
  
		- **3_upd_base_next_pos** : ベースの移動位置を決定する関数
  
		- **4_upd_kinematic_feasibility** : 次の遊脚・把持点・ベース位置が逆運動学的に可能か判定する関数
  
		- **5_upd_gait_history : 歩容計画の結果を保存する関数
  
		- `ini_gait.m` : 歩容計画に関わる変数の初期設定をおこなう関数
  
		- `upd_gait_planning.m` : 歩容計画の関数群のおおもとであり，`main_sim`のループ内で実行される歩容計画パラメータ更新の関数
  
	- **global path plan** : 大局的経路計画に関する関数群を収容．
  
	- **local_path_plan** : 局所的経路計画に関する関数を収容
  
	- **motion** : ベースの動き・脚の軌道を生成・制御する関数群を収容．
  
	- **tools** : 順運動学の関数と，歩容に関わるmap情報を取得する関数が収容
  
	- `upd_grasp_detach.m` : グリッパの滑落を判定して更新する関数
  
- **engine** : ロボットにかかる外力の計算や外部との接触モデルに関わる関数．

- **environment**

	- **climbing_holds_map** : ボルダリングホールドを用いた地形情報の保管場所

	- **grid_map** : グリッド形式の生成した地形情報の保管場所
  
	- **terrain_map** : 実地形マップの保管場所
  
	- **uneven_map** : フラクタルを利用した不整地マップの保管場所
  
	- `map_grid_desingner.m` : グリッド形式のマップの生成をする実行ファイル
  
	- 地形情報・シミュレーション環境に関するinitialize関数群 : `ini_environment.m`, `ini_surface.m`, `ini_graspable_points.m`
  
	- `map_format_conversion.m` : RealSenseカメラで取得してPCDから3×N行列とした点群情報を，座標変換・線形補間してグリッド形式に変換する関数．

- **misc** : 未分類の関数群．データ保存に関わるsave関数などが置かれている．

- **robot** : ロボットモデルの定義やその初期設定・解析を行う関数を収容．

	- **LP_analysis** : ロボットのリンクパラメータ(LP)から脚先可動範囲を描画する関数
  
	- `ini_12DOF_SV.m`, `ini_nDOF_SV.m`: 設定されたロボットの自由度にもとづき状態変数(SV)を設定する関数
	
	- `ini_(Robot名)関数群` : 各ロボットのLPを内部で指定してメインに読み込む
	
	- `ini_joint_angle.m`
	
	- `ini_reachable_area.m`
	
	- `ini_robot.m`

- **sensing** : カメラやセンサにもとづくセンシングに関連するini_,upd_関数群を主要

- **visualization** : 全てのvisualize関数群の保存場所

- `ini_path.m` : `main_sim.m`の実行のために，実行しているフォルダの現在位置確認と，パス設定をおこなう．

- ``main_sim.m` : ClimbLabの実行ファイル

<br>

以上がClimbLabのフォルダ構成の概要である．
