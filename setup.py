from setuptools import setup

setup(
      name="mbuild_CLP",
      install_requires="mbuild",
      entry_points={
                    "mbuild.plugins":[ "CLP = mbuild_CLP.mbuild_CLP:CLP_box"
                        ]
                    },
                    py_modules=["mbuild_CLP"],
                        )
