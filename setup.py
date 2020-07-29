from setuptools import setup

setup(
<<<<<<< HEAD
      name="mbuild_CLP",
      install_requires="mbuild",
      entry_points={
                    "mbuild.plugins":[ "CLP = mbuild_CLP.mbuild_CLP:CLP_box"
                        ]
                    },
                    py_modules=["mbuild_CLP"],
                        )
=======
    name="mbuild_ona",
    install_requires="mbuild",
    entry_points={"mbuild.plugins": ["DNA = mbuild_ona.mbuild_ona:DNA"]},
    py_modules=["mbuild_ona"],
)
>>>>>>> 704e82ccbfc1088140a9cd18d7b6c468559b3963
