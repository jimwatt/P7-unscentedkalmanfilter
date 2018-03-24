# Sensor Fusion with the Unscented Kalman Filter

This project is in conjunction with the Udacity Self-Driving Car course.  In particular, the goal of this project is to use the Unscented Kalman Filter to fuse simulated radar and laser measurements and estimate the position of a cyclist.

* More detail about the approach is provided in [writeup.md](./writeup.md).

## Getting Started

### Prerequisites

```
CMake
uWS		(Socket for interfacing with Udacity simulator)	
```

### Installing

1. Clone this project from GitHub:

```
git clone https://github.com/jimwatt/P7-unscentedkalmanfilter.git
```

2. Change into the `build` directory, and compile with `CMake`.

   ```
   cd build
   cmake ..
   make
   ```

## Running the Code

* Run the executable from the `build` directory:
  ```./UnscentedKF```
* The executable will wait until measurements are provided by the simulator through the uWS socket.





## Authors

* **James Watt**

<!--## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details-->

## Acknowledgments
This project is a submission to the Udacity Self-Driving Car nanodegree:

* ```
  https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project.git
  ```

