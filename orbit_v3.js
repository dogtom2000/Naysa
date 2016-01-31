var planet = ["Earth", 3.986e14, 6.371e6, 1, 7, 1.40e5, 28.97, 24];

//Planet Name
//Gravitaional parameter (m3/s2)
//Radius (m)
//Pressure (atm)
//Atmospheric Scale (dimensionless)
//Atmospheric Height (m)
//Atmospheric Molecular Weight (g/mol)
//Rotational period (hours)

var rocket = ["Atlas_MkII", [900, 15000, 250000, 300, 0.2, 5], [300, 5000, 65000, 400, 0.2, 2.5]];

// Rocket Name
// Stag Properties:
// [Dry Mass (kg)
// Fuel Mass (kg)
// Max Thrust (N)
// Specific Impluse (s)
// Drag Coefficient (dimensionless)
// Diamter (m)]

//Execute from here
var output = orbitPlanet(planet,rocket,180,Math.PI / 2);

function orbitPlanet(planet, rocket, orbit, phi){
    //construct planet and rocket objects
    var Planet = new planetConst(planet);
    var Rocket = new rocketConst(rocket);
    
    //define dynamic rocket variables
    var time = [];
    var acceleration = [[0,0,0]];
    var velocity = [[0,0,0]];
    var position = [[0,0,0]];

    //define rocket variables
    var rocketArea = 0;
    var maxFuel = 0;
    var rocketStage = [];
    var stageMass = 0;
    var stage = 0;

    //construct rocket stage arrays
    for (var i = 0; i < Rocket.massDry.length; i++){
        
        stageMass = 0;
        for (var j = i; j < Rocket.massDry.length; j++){
            stageMass += Rocket.massDry[j] + Rocket.massFuel[j]
        }
        rocketStage[i] = [stageMass, Rocket.massFuel[i], Rocket.thrust[i], Rocket.Isp[i], Rocket.cD[i], Rocket.diameter[i]];
        
    }
    
    //calculate initial velocity
    velocity[0] = [0, Planet.radius * Math.sin(phi) * Planet.rotation, 0];

    //calculate initial position
    position[0] = [Planet.radius, 0, phi];
    var timeStep = 0;
    for (var k = 0; k < 2; k++){
    var output = burn(rocketStage, stage, acceleration, velocity, position, orbit, Planet, timeStep);
    
    var apoFlag = output[0];
    var stopFlag = output[1];;
    stage = output[2];;
    rocketStage[stage][2] = output[3];;
    acceleration.push(output[4]);
    velocity.push(output[5]);
    position.push(output[6]);
    timeStep++;
    console.log(output[4])
    }

}

function burn(rocketStage, stage, acceleration, velocity, position, orbit, Planet, timeStep){
    //define flow control variables
    var apoFlag = 0;
    var stopFlag = 0;
    
    if (rocketStage[stage][1] > 0){
        var currentStage = rocketStage[stage];
    } else if(rocketStage.length > stage + 1) {
        stage++;
        var currentStage = rocketStage[stage]
    } else {
        var currentStage = 0;
    }
    
    var deltaTime = 10;
    
    var stageTotalMass = currentStage[0];
    var stageFuelMass = currentStage [1];
    var stageDryMass = stageTotalMass - stageFuelMass;
    var stageThrust = currentStage[2];
    var stageIsp = currentStage[3];
    
    
    
    //calculate attributes
    var surfaceVelocity = [velocity[timeStep][0], velocity[timeStep][1] - position[timeStep][0] * Math.sin(position[timeStep][2]) * Planet.rotation, velocity[timeStep][2]];
    var surfaceSpeed = magni(surfaceVelocity);
    var orbitalSpeed = magni(velocity[timeStep]);
    var angMom = cross(polarToCartesian(position[timeStep]), velocity[timeStep]);
    var vCrossH = cross(velocity[timeStep], angMom);
    var ecc = [vCrossH[0] / Planet.gravPara - 1,vCrossH[1] / Planet.gravPara, vCrossH[2] / Planet.gravPara];
    var orbitEnergy = Math.pow(orbitalSpeed, 2) / 2 - Planet.gravPara / magni(polarToCartesian(position[timeStep]));
    var smAxis = - Planet.gravPara / 2 / orbitEnergy;
    var apoa = smAxis * (1 + Math.abs(magni(ecc)));
    var peri = smAxis * (1 - Math.abs(magni(ecc)));
    
    //centrifugal and gravitational accelerations
    var cent = (Math.pow(velocity[timeStep][1],2) + Math.pow(velocity[timeStep][2],2)) / position[timeStep][0];
    var grav = Planet.gravPara / Math.pow(position[timeStep][0],2);
    
    //drag calculations
    if (position[timeStep][0] > Planet.atmHeight + Planet.radius || surfaceSpeed === 0){
        drag = [0,0];
    } else {
        var density = 0.042295 * Planet.pressure * Math.exp(-(position[timeStep][0] - Planet.radius) / Planet.atmScale / 1000) * Planet.atmWeight;
        var dragMagnitude = currentStage[4] * density * Math.pow(surfaceSpeed, 2) * Math.PI / 4 * Math.pow(currentStage[5],2) / 2 / stageTotalMass;
        var drag = surfaceVelocity.map(function(num) {return num * dragMagnitude / surfaceSpeed});
    }
    
    if (apoFlag === 0){
        var pitch = Math.PI / 2;
    } else if (velocity[1] > 0) {
        pitch = Math.asin(drag[0] / stageThrust * stageTotalMass);
    } else {
        pitch = Math.asin((drag[0] + grav) / stageThrust * stageTotalMass);
    }
    
    var thrust = [stageThrust * Math.sin(pitch) / stageTotalMass, stageThrust * Math.cos(pitch) / stageTotalMass, 0];
    
    acceleration[timeStep] = [thrust[0] + drag[0] + cent - grav, thrust[1] - drag[1], 0];
    
    for (var i = 0; i < 3; i++){
        velocity[timeStep][i] += acceleration[timeStep][i] * deltaTime;
    }
    

    
    if (timeStep === 0){
        position[timeStep][0] += velocity[timeStep][0] / 2;
        position[timeStep][1] += velocity[timeStep][1] / position[timeStep][0] / 2;
        position[timeStep][2] += velocity[timeStep][2] / position[timeStep][0] / 2;                  
    } else {
        position[timeStep][0] += (velocity[timeStep][0] + velocity[timeStep-1][0]) / 2;
        position[timeStep][1] += (velocity[timeStep][1] + velocity[timeStep-1][1]) / position[0] / 2;
        position[timeStep][2] += (velocity[timeStep][2] + velocity[timeStep-1][2]) / position[0] / 2;  
    }
    
    

    stageFuelMass -= stageThrust / stageIsp / 9.81 * deltaTime;
    
    if (apoa > orbit + Planet.radius && apoFlag === 0){
            
        apoFlag = 1;

    }
    
    if (cent > grav){
        return stopFlag = 1;
    }
    

    return [apoFlag, stopFlag, stage, stageFuelMass, acceleration, velocity, position];
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
};

function rocketConst(rocketArray) {
    this.name = rocketArray[0];
    var massDry = [];
    var massFuel = [];
    var thrust = [];
    var Isp = [];
    var cD = [];
    var diameter = [];
    
    for (var i = 1; i < rocketArray.length; i++){
    massDry.push(rocketArray[i][0]);
    massFuel.push(rocketArray[i][1]);
    thrust.push(rocketArray[i][2]);
    Isp.push(rocketArray[i][3]);
    cD.push(rocketArray[i][4]);
    diameter.push(rocketArray[i][5]);
    }
    
    this.massDry = massDry;
    this.massFuel = massFuel;
    this.thrust = thrust;
    this.Isp = Isp;
    this.cD = cD;
    this.diameter = diameter;

    this.position = [0,0,0];
    this.velocity = [0,0,0];
};

function magni(vector){
    
    return Math.pow(Math.pow(vector[0],2) + Math.pow(vector[1],2) + Math.pow(vector[2],2),0.5);
   
}

function polarToCartesian(vector){
    
    var x = vector[0] * Math.cos(vector[1]) * Math.sin(vector[2]);
    var y = vector[0] * Math.sin(vector[1]) * Math.sin(vector[2]);
    var z = vector[0] * Math.cos(vector[2]);
    
    //console.log(vector, [x, y, z])
    return [x, y, z];
    
}

function cross(u, v){
    
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    
}
