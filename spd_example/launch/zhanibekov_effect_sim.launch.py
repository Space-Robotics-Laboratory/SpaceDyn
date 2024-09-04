import os 

from ament_index_python.packages import get_package_share_directory

from launch import LaunchDescription
from launch_ros.actions import Node

import xacro

pkg_dir = get_package_share_directory('spd_example')
launch_file_path = os.path.join(pkg_dir, 'launch', 'zhanibekov_effect_sim.launch.py')
xacro_path = os.path.join(pkg_dir, 'urdf', 'zhanibekov_effect_model.urdf.xacro')
urdf_path = os.path.join(pkg_dir, 'urdf', 'zhanibekov_effect_model.urdf')
rviz_path = os.path.join(pkg_dir, 'config', 'zhanibekov_effect_sim.rviz')

 # make urdf from xacro
# load xacro
doc = xacro.process_file(xacro_path)
# make urdf
robot_desc = doc.toprettyxml(indent=' ')
# export urdf to urdf path
f = open(urdf_path, 'w')
f.write(robot_desc)
f.close()

def generate_launch_description():
  rviz2 = Node(
    package="rviz2",
    executable="rviz2",
    name="rviz2",
    arguments=["-d", rviz_path]
  )

  zhanibekov_effect_sim = Node(
    package="spd_example",
    executable="zhanibekov_effect_sim",
    name="zhanibekov_effect_sim"
  )

  robot_state_publisher_node = Node(
    package='robot_state_publisher',
    executable='robot_state_publisher',
    name='robot_state_publisher',
    arguments=[urdf_path]
  )

  return LaunchDescription([
    rviz2,
    zhanibekov_effect_sim,
    robot_state_publisher_node
  ])
