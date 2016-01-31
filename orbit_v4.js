function orbitEarth(planet, rocket){
//[things to add
//trajectory function with different trajectory options
//inclination calculation
//phi options (might be buggy)]

//var planet = ["Earth", 3.986e14, 6.371e6, 1, 7, 1.40e5, 28.97, 24];

//[Planet Name
//Gravitaional parameter (m3/s2)
//Radius (m)
//Pressure (atm)
//Atmospheric Scale (dimensionless)
//Atmospheric Height (m)
//Atmospheric Molecular Weight (g/mol)
//Rotational period (hours)]

//var rocket = ["Saturn V", [130000, 2160000, 2160000, 34000000, 263, 0.2, 10.1], [40100, 456100, 456100, 4400000, 412, 0.2, 10.1], [74300, 109500, 109500, 1000000, 412, 0.2, 6.6]];
//var rocket = ["Atlas_MkII", [1500, 15000, 15000, 320000, 250, 0.2, 5], [500, 5000, 5000, 150000, 300, 0.2, 2.5], [500, 3000, 3000, 65000, 300, 0.2, 2.5]];
// [Rocket Name
// Stag Properties:
// Dry Mass (kg)
// Fuel Mass (kg)
// Max Thrust (N)
// Specific Impluse (s)
// Drag Coefficient (dimensionless)
// Diamter (m)]

//Execute from here
var output = orbitPlanet(planet, rocket, 173000, 0, Math.PI / 2);
//console.log(output[0])
return output;
//console.log(output);

function orbitPlanet(planet, rocket, orbit, theta, phi){
    
    //construct planet and rocket objects
    var Planet = new planetConst(planet);
    var Rocket = new rocketConst(rocket);
    
    var Mass1 = [];
    var dV1 = 0;
    for (var i = 0; i < Rocket.massDry.length; i++){
    Mass1[i] = 0;
        for (var j = i; j < Rocket.massDry.length; j++){
            Mass1[i] += Rocket.massDry[j] + Rocket.massFuel[j]
        }
        dV1 += Math.log(Mass1[i] / (Mass1[i] - Rocket.massFuel[i])) * Rocket.Isp[i] * 9.81
    }
    //console.log(dV1)
    //define flow control variables
    var t = 0;
    var deltaTime = 1;
    var maxRun = 1500;
    var stopFlag = 0;
    var apoFlag = 0;
    var outputMessage = "";
    
    //construct stage definitions
    var rocketStage = stageConst(Rocket);
    
    //define rocket variables
    var stage = 0;
    var time = [0];
    var heading = []; 
    var acceleration = [[0,0,0]]; //[d2r/dr2, d2theta/dt2, d2phi/dt2]
    var velocity = [[0,0,0]]; //[dr/dt, dtheta/dt, dphi/dt]
    var position = [[0,0,0]]; //[r, theta, phi]
    var positionAddLast = [0, Math.sin(phi) * Planet.rotation, 0];
    
    //intial values
    velocity[0] = [0, Planet.radius * Math.sin(phi) * Planet.rotation, 0];
    position[0] = [Planet.radius, theta, phi];

    var output = burn();
    //console.log(outputMessage);
    //console.log(output)
    return output;
    
    function burn(){
        for (var run = 0; run < maxRun + 1; run++){
            if (run == maxRun){
                outputMessage = "max loops reached";
                stopFlag = 1;
            }
            
            //define current stage or break if out of fuel
            if (rocketStage[stage][1] > 0){
                var currentStage = rocketStage[stage];
            } else if(rocketStage.length > stage + 1) {
                rocketStage[stage][1] = 0;
                stage++;
                currentStage = rocketStage[stage];
            } else {
                outputMessage = "out of fuel";
                stopFlag = 1;
            }

            //define stage variables
            var stageMassTotal = currentStage[0]; //does not change
            var stageMassFuel = currentStage [1]; //changes on each time step
            var stageMassDry = stageMassTotal - currentStage[2]; //does not change
            var stageMassCurrent = stageMassDry + stageMassFuel; //changes on each time step
            var stageThrust = currentStage[3];
            var stageIsp = currentStage[4];
            var stageBurnRate = stageThrust / stageIsp / 9.81;

            if (stageThrust / (stageMassCurrent * 9.81) < 1 && stage == 0){
                outputMessage = "first stage TWR must be greater than 1";
                stopFlag = 1;
            }
        
            //calculate orbital variables
            var surfaceVelocity = arrayAdd(velocity[t], [0, - Planet.radius * Math.sin(phi) * Planet.rotation, 0]);
            var surfaceSpeed = magn(surfaceVelocity);
            var speed = magn(velocity[t]);
            var angularMomentum = cross([position[t][0], 0, 0], velocity[t]);
            var eccVector = arraySub(arrayDiv(cross(velocity[t], angularMomentum), Planet.gravPara), [1, 0, 0]);
            var ecc = magn(eccVector);
            var orbitalEnergy = Math.pow(speed, 2) / 2 - Planet.gravPara / position[t][0];
            var semiMajorAxis = -Planet.gravPara / 2 / orbitalEnergy;
            var apoapsis = semiMajorAxis * (1 + Math.abs(ecc));
            var periapsis = semiMajorAxis * (1 - Math.abs(ecc));

            //calculate acceleration [d2r/dr2, d2theta/dt2, d2phi/dt2]
            var centripetalAcceleration = [(Math.pow(velocity[t][1],2) + Math.pow(velocity[t][2],2)) / position[t][0], 0, 0];
            var gravityAcceleration = [-Planet.gravPara / Math.pow(position[t][0],2), 0, 0];
            if (surfaceSpeed === 0 || position[t][0] > Planet.atmHeight + Planet.radius){
                var dragAcceleration = [0,0,0];
            } else {
                var density = 0.042295 * Planet.pressure * Math.exp(-(position[t][0] - Planet.radius) / Planet.atmScale / 1000) * Planet.atmWeight;
                var dragFraction = -currentStage[5] * density * Math.pow(surfaceSpeed, 2) * currentStage[6] / 2 / stageMassCurrent / surfaceSpeed;
                dragAcceleration = arrayMul(surfaceVelocity, dragFraction);
            }
            
            //heading calculations: r of position vector is alligned with x of heading vector 
            //[x, y, z] = [dr/dt, dtheta/dt, dphi/dt] = [rcos(theta)sin(phi), rsin(theta)sin(phi), rcos(phi)]
            //theta = 0, phi = pi/2 = pointed in r direction
            //heading = [theta, phi]
            
            if (position[t][0] - Planet.radius < Planet.atmHeight / 10 + 250){
                heading = [0, Math.PI / 2];
            } else {
                heading = [Math.PI / 2 * (position[t][0] - Planet.radius - Planet.atmHeight / 10 + 250) / (orbit - Planet.atmHeight / 10 + 250), Math.PI / 2];
            }
            
            var thrustVector = [Math.cos(heading[0]) * Math.sin(heading[1]), Math.sin(heading[0]) * Math.sin(heading[1]), Math.cos(heading[1])];
            var thrustAcceleration = arrayMul(thrustVector, stageThrust / stageMassCurrent);
          
            t++;
            
            acceleration[t] = arrayAddPlus(centripetalAcceleration, gravityAcceleration, dragAcceleration, thrustAcceleration);
            
            apoFlag = 0;
            if (apoapsis - Planet.radius > orbit && position[t - 1][0] - Planet.radius < Planet.atmHeight){

                stageBurnRate = magn(dragAcceleration) * stageMassCurrent / stageIsp / 9.81;
                acceleration[t] = arrayAdd(centripetalAcceleration, gravityAcceleration);
                apoFlag = 1;
                
            } else if (apoapsis - Planet.radius > orbit && position[t - 1][0] - Planet.radius >= Planet.atmHeight) {
                
                acceleration[t] = arrayAdd(centripetalAcceleration, gravityAcceleration);
                apoFlag = 2;
                
            }

            var accelTime = 20 / magn(acceleration[t]);
            
            var fuelCount = 10;
            
            while (fuelCount < stageMassFuel){
                fuelCount *= 10;
            }
            
            var burnCount = 10;
            
            while (burnCount < accelTime * stageBurnRate){
                burnCount *= 10;
            }
            
            var burnTick = burnCount / 10;
            
            if (burnTick > fuelCount / 10){
                burnTick = fuelCount / 10;
            }

            deltaTime = burnTick / stageBurnRate;
            
            if (stageBurnRate < 1 && apoFlag == 1){
                stageBurnRate = 1;
                deltaTime = 1 / stageBurnRate;
            }
            time[t] = time[t - 1] + deltaTime;
            
            velocity[t] = arrayAdd(velocity[t - 1], arrayMul(acceleration[t], deltaTime));
            
            var positionAdd = [velocity[t][0], velocity[t][1] / position[t - 1][0], velocity[t][2] / position[t - 1][0]];
            
            var positionAddAve = arrayMul(arrayAdd(positionAdd, positionAddLast), 0.5);
            
            position[t] = arrayAdd(position[t - 1], arrayMul(positionAddAve, deltaTime)); 

            positionAddLast = positionAdd;
            
            if (position[t][0] < planet.Radius){
                outputMessage = "you have crashed, everyone is dead";
                stopFlag = 1;
            }

            if (apoFlag !== 2){
            currentStage[1] = Math.round(currentStage[1] - stageBurnRate * deltaTime);
            if (currentStage[1] < 0){
                currentStage[1] = 0;
            }
            
            
            }
            
            if (position[t][0] > orbit + Planet.radius){
                
                var dvRequired = Math.pow(Planet.gravPara / (orbit + Planet.radius), 0.5) - magn(velocity[t]);
                var stageDv = Math.log(stageMassCurrent / stageMassDry) * stageIsp * 9.81;

                while (dvRequired > 0 && rocketStage.length > stage){
                    
                    if (stageThrust / (stageMassCurrent * 9.81) < 0.5){
                        outputMessage = "not enough thrust to complete burn";
                        stopFlag = 1;
                        break;
                    }
                    
                    if (dvRequired > stageDv){
                        dvRequired -= stageDv;
                        stageDv = 0;
                        currentStage[1] = 0;
                        rocketStage[stage][1] = currentStage[1];
                    } else {
                        stageDv -= dvRequired;
                        dvRequired = 0;
                        currentStage[1] = Math.round(Math.exp(stageDv / stageIsp / 9.81) * stageMassDry - stageMassDry);
                        rocketStage[stage][1] = currentStage[1];
                    }
                    stage++;
                    if (dvRequired > 0 && rocketStage.length > stage){
                        currentStage = rocketStage[stage];
                        stageMassTotal = currentStage[0]; //does not change
                        stageMassFuel = currentStage [1]; //changes on each time step
                        stageMassDry = stageMassTotal - currentStage[2]; //does not change
                        stageMassCurrent = stageMassDry + stageMassFuel; //changes on each time step
                        stageThrust = currentStage[3];
                        stageIsp = currentStage[4];
                        stageDv = Math.log(stageMassCurrent / stageMassDry) * stageIsp * 9.81;
                    }
                }
                
                if (stopFlag !== 1){
                    
                    outputMessage = "dV required for orbit: " + dvRequired + ", dV remaining: " + stageDv + " on stage" + stage;
                    stopFlag = 1;
                }

            }
            
            if (stopFlag === 1){
                rocket = stageRocketConst(Rocket, rocketStage);
                //return [rocket, velocity[t], position[t], run];
                return [time ,position, velocity, acceleration, rocket, run, outputMessage];
            }

        } 
    }
}

function planetConst(planetArray) {
    this.name = planetArray[0];
    this.gravPara = planetArray[1];
    this.radius = planetArray[2];
    this.pressure = planetArray[3];
    this.atmScale = planetArray[4];
    this.atmHeight = planetArray[5];
    this.atmWeight = planetArray[6];
    this.rotation = 2 * Math.PI / (planetArray[7] * 3600);
}

function rocketConst(rocketArray) {
    this.name = rocketArray[0];
    var massDry = [];
    var massFuel = [];
    var massFuelCapacity = [];
    var thrust = [];
    var Isp = [];
    var cD = [];
    var diameter = [];
    
    for (var i = 1; i < rocketArray.length; i++){
    massDry.push(rocketArray[i][0]);
    massFuel.push(rocketArray[i][1]);
    massFuelCapacity.push(rocketArray[i][2]);
    thrust.push(rocketArray[i][3]);
    Isp.push(rocketArray[i][4]);
    cD.push(rocketArray[i][5]);
    diameter.push(rocketArray[i][6]);
    }
    
    this.massDry = massDry;
    this.massFuel = massFuel;
    this.massFuelCapacity = massFuelCapacity;
    this.thrust = thrust;
    this.Isp = Isp;
    this.cD = cD;
    this.diameter = diameter;

}

function stageConst(Rocket){
    
    var rocketStage = [];
    for (var i = 0; i < Rocket.massDry.length; i++){
    var stageMass = 0;
        for (var j = i; j < Rocket.massDry.length; j++){
            stageMass += Rocket.massDry[j] + Rocket.massFuel[j]
        }
        rocketStage[i] = [stageMass, Rocket.massFuel[i], Rocket.massFuel[i], Rocket.thrust[i], Rocket.Isp[i], Rocket.cD[i], Math.PI / 4 *Math.pow(Rocket.diameter[i],2)];
    }
    return rocketStage;
    
}

function stageRocketConst(Rocket, rocketStage){
    
    var rocket = [Rocket.name];
    for (var i = 0; i < rocketStage.length; i++){
        rocket.push([Rocket.massDry[i], rocketStage[i][1], Rocket.massFuelCapacity[i], Rocket.thrust[i], Rocket.Isp[i], Rocket.cD[i], Rocket.diameter[i]]);
    }
    return rocket;
}

function arrayAdd(array, val){
    
    var arrayNew = [];
    if (Array.isArray(val)){
        if (array.length !== val.length){
            return "array length mismatch"
        } else {
            for (var i = 0; i < array.length; i++){
                arrayNew[i] = array[i] + val[i];
            }            
        } 
    } else {
        for (var i = 0; i < array.length; i++){
            arrayNew[i] = array[i] + val;
        }
    }
    return arrayNew;

}

function arraySub(array, val){
    
    var arrayNew = [];
    if (Array.isArray(val)){
        if (array.length !== val.length){
            return "array length mismatch"
        } else {
            for (var i = 0; i < array.length; i++){
                arrayNew[i] = array[i] - val[i];
            }            
        } 
    } else {
        for (var i = 0; i < array.length; i++){
            arrayNew[i] = array[i] - val;
        }
    }
    return arrayNew;

}

function arrayMul(array, val){
    
    var arrayNew = [];
    if (Array.isArray(val)){
        if (array.length !== val.length){
            return "array length mismatch"
        } else {
            for (var i = 0; i < array.length; i++){
                arrayNew[i] = array[i] * val[i];
            }            
        } 
    } else {
        for (var i = 0; i < array.length; i++){
            arrayNew[i] = array[i] * val;
        }
    }
    return arrayNew;

}

function arrayDiv(array, val){
    
    var arrayNew = [];
    if (Array.isArray(val)){
        if (array.length !== val.length){
            return "array length mismatch"
        } else {
            for (var i = 0; i < array.length; i++){
                arrayNew[i] = array[i] / val[i];
            }            
        } 
    } else {
        for (var i = 0; i < array.length; i++){
            arrayNew[i] = array[i] / val;
        }
    }
    return arrayNew;

}

function arrayAddPlus(){
    
    var newArray = [];
    for (var i = 0; i < arguments[0].length; i++){
        newArray[i] = 0;
        for (var j = 0; j < arguments.length; j++){
           newArray[i] += arguments[j][i]; 
        }
    }
  return newArray;
  
}

function magn(vector){
    
    if (vector.length == 2){
        return Math.pow(Math.pow(vector[0],2) + Math.pow(vector[1],2),0.5);
    } else if (vector.length == 3){
        return Math.pow(Math.pow(vector[0],2) + Math.pow(vector[1],2) + Math.pow(vector[2],2),0.5);
    }
   
}

function s2c(vector){
    
    var x = vector[0] * Math.cos(vector[1]) * Math.sin(vector[2]);
    var y = vector[0] * Math.sin(vector[1]) * Math.sin(vector[2]);
    var z = vector[0] * Math.cos(vector[2]);
    return [x, y, z];
    
}

function c2s(vector){
    
    var r = Math.pow(Math.pow(vector[0],2) + Math.pow(vector[1],2) + Math.pow(vector[2],2),0.5);
    var t = Math.atan(vector[1] / vector[0]);
    var p = Math.acos(vector[2] / r);
    return [r, t, p];
}

function cross(u, v){
    
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    
}
}