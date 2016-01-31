var planet = ["Earth", 398600, 6371, 1, 7, 140, 28.97, 465];

//Planet Name
//Gravitaional parameter (km3/s2)
//Radius (km)
//Pressure (atm)
//Atmospheric Scale (dimensionless)
//Atmospheric Height (km)
//Atmospheric Molecular Weight (g/mol)
//Tangential Velocity at Equator (m/s)

var rocket = ["Atlas_MkII", 1000, 23000, 270000, 600, 0.2, 5];

// Rocket Name
// Dry Mass (kg)
// Fuel Mass (kg)
// Max Thrust (N)
// Specific Impluse (s)
// Drag Coefficient (dimensionless)
// Diamter (m)

//Execute from here
for (var j = 1; j < 5; j++)  {
    
    var output = burn(planet,rocket,180);
    
    rocket[4] -= 100;
}


function burn(planet, rocket, orbit){
    //construct planet and rocket objects
    var Planet = new planetConst(planet);
    var Rocket = new rocketConst(rocket);
    
    //define flow control variables
    var apoFlag = 0;
    var stopFlag = 0;
    
    //define universal constants
    var gravConst = 9.81; //m/s

    //define rocket specific consants
    var rocketArea = Math.PI / 4 * Math.pow(Rocket.diameter,2); //m2
    var maxFuel = Rocket.massFuel; //kg
    
    //define rocket variables
    var velocity = [0,0];
    var position = [0,0];
    var apoapsis = 0;
    var periapsis = 0;
    var eccentricity = 1;
    var smAxis = 0;
    var velocityLast = [0,0];

    for (var i = 1; i < 15000; i++){
        
        if (stopFlag === 0){
            accelerationCalc(maxFuel / (Rocket.thrust / Rocket.Isp / gravConst) / 10000);
        } else {
            console.log(Rocket)
            break;
        }
    }

    function accelerationCalc(deltaTime){
        //Calculations
        var massShip = Rocket.massDry + Rocket.massFuel;
        var shipRadius = position[1] + Planet.radius;
        var speed = vectorMag(velocity);
        var orbitalVelocity = [velocity[0] + Planet.velEquator, velocity[1]];
        var orbitalSpeed = vectorMag(orbitalVelocity);
        eccentricity = (Math.pow(orbitalVelocity[0],2) * shipRadius + orbitalVelocity[0] * orbitalVelocity[1] * shipRadius) / Planet.gravPara / 1e6 - 1;
        smAxis = -Planet.gravPara / 2 / (Math.pow(orbitalSpeed, 2) / 2e6 - Planet.gravPara / shipRadius);
        apoapsis = smAxis * ( 1 + Math.abs(eccentricity));
        periapsis = smAxis * ( 1 - Math.abs(eccentricity));

        //centrifugal and gravitational accelerations
        var cent = Math.pow(orbitalVelocity[0],2) / (Planet.radius + position[1]) / 1000;
        var grav = Planet.gravPara / Math.pow((Planet.radius + position[1]),2) * 1000;
        
        //drag calculations
        if (position[1] > Planet.atmHeight || speed === 0){
            drag = [0,0];
        } else {
            var density = 0.042295 * Planet.pressure * Math.exp(-position[1] / Planet.atmScale) * Planet.atmWeight;
            var dragMagnitude = Rocket.cD * density * Math.pow(speed, 2) * rocketArea / 2 / massShip;
            var drag = velocity.map(function(num) {return num * dragMagnitude / speed});
        }
        
        if (apoFlag === 0){
            var pitch = Math.PI / 2;
        } else if (velocity[1] > 0) {
            pitch = Math.asin(drag[1] / Rocket.thrust * massShip);
        } else {
            pitch = Math.asin((drag[1] + grav) / Rocket.thrust * massShip);
        }
    
        var thrust = [Rocket.thrust * Math.cos(pitch) / massShip, Rocket.thrust * Math.sin(pitch) / massShip];

        velocity[0] += (thrust[0] - drag[0]) * deltaTime;
        velocity[1] += (thrust[1] - drag[1] + cent - grav) * deltaTime;
        position[1] += (velocityLast[1] + velocity[1]) / 2 * deltaTime / 1000;
        
        velocityLast = velocity;
        
        Rocket.massFuel -= Rocket.thrust / Rocket.Isp / gravConst * deltaTime;
        
        if (apoapsis > orbit + Planet.radius && apoFlag === 0){
            
            apoFlag = 1;

        }
        
        if (cent > grav){
            console.log("orbit achieved")
            Rocket.position = [0,orbit];
            Rocket.velocity = [1000 * Math.pow(Planet.gravPara / (Planet.radius + orbit),0.5),0];
            return stopFlag = 1;
        }
        
        if (Rocket.massFuel <= 0){
            console.log("out of fuel")
            console.log("apoapsis", (apoapsis), "periapsis = ", (periapsis))
            Rocket.position = [0,apoapsis - 6371];
            Rocket.velocity = velocity;
            return stopFlag = 1;
            
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
    this.velEquator = planetArray[7];
};

function rocketConst(rocketArray) {
    this.name = rocketArray[0];
    this.massDry = rocketArray[1];
    this.massFuel = rocketArray[2];
    this.thrust = rocketArray[3];
    this.Isp = rocketArray[4];
    this.cD = rocketArray[5];
    this.diameter = rocketArray[6]
    this.position = [0,0];
    this.velocity = [0,0]
};

function vectorMag(vector){
    
    var magnitude = Math.pow((Math.pow(vector[0],2) + Math.pow(vector[1],2)),0.5);
    return magnitude;
    
}
