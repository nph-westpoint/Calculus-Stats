# conoculator2.0
Capturing Calculus Measurements from Discrete Data

### General Information
This repository is written for the nutrition precision health researchers who are currently using discrete data to study the Mixed Meal Response Curves on homeostasis of blood sugar levels in the body. The theory is that individuals with post-prandial response to meals that deviate more from homeostasis than normal have a risk factor.  
Glucose values are measured after a pre-determined amount of fasting has occurred in order to establish the normal homeostasis value for an individual. A standard meal is used to test the body's response to glucose. The blood glucose levels are recorded at specified times and recorded as (time, glucose level) points that can be used with the methods located in the repository to gain insight about the health of the patient.

### Current Idea - Incremental Area Under the Curve (iAUC) 
Currently, iAUC is the only calculation used to determine the health of the patient. iAUC is the difference between the total area under the curve and normal homeostasis value.
![image](https://github.com/nph-westpoint/conoculator/assets/142028542/d4e21d31-8022-4832-b030-dbafb3e6802c)

Note that the blue shade is not the total area, but the deviation from the time that the meal was given which is why it is referred to as the incremental area under the curve. Also, a standard amount of time to measure is 2 hours which you can see from the picture.

### Purpose
The authors of the conoculator propose to use more calculus measures than iAUC as a means of characterizing the shape of a post-prandial glucose curve. Additionally, we want to make the calculation of these calculus based statistics easy for anyone with discrete data regardless of the interval that they used while also allowing them the ability to impute (fill in missing data) time series data without getting technical in the mathematics.



### Data Dictionary for Results
- auc180: The area under the curve for 180 minutes.
- auc240: The area under the curve for 240 minutes.
- iauc180: The incremental area under the curve for 180 minutes. Note that if your data stops before 180 minutes, that this measure will be erroneous.
- iauc240: The incremental area under the curve for 240 minutes. Note that if your data stops before 240 minutes, that this measure will be erroneous.
- max_inc_time:
- max_inc_val:
- max_dec_time:
- max_dec_val:
- abs_max_time:
- abs_max_val:
- inf_pt_time:
- inf_pt_val:
- tot_crit_pts:
- local_max:
- local_min:


### Documentation
All of the numerical analysis techniques (cubic splines and 5-point derivatives) used in this software was developed using the Numerical Analysis book by Burden and Faires. The paper used to base many calculations on is the Zeevi paper.
- Burden, Richard and J. Douglas Faires Numerical Analysis (9th ed.). (2011). Brooks/Cole.
- Personalized Nutrition by Prediction of Glycemic Responses. Zeevi D, Korem T, Zmora N, Israeli D, Rothschild D, Weinberger A, Ben-Yacov O, Lador D, Avnit-Sagi T, Lotan-Pompan M, Suez J, Mahdi JA, Matot E, Malka G, Kosower N, Rein M, Zilberman-Schapira G, Dohnalová L, Pevsner-Fischer M, Bikovsky R, Halpern Z, Elinav E, Segal E. Cell. 2015 Nov 19;163(5):1079-1094. doi: 10.1016/j.cell.2015.11.001.
