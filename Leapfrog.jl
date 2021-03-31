"""
Function for main numerical method

    Takes input array of all vortecies and applies numerical methods on them
    By modifying parameters, the output may be calculated.

    Requires use of "ConfParser" and "Plots" libraries.
"""


#data type for making a vortex with all associated parameters and variables
mutable struct Vortex
        position #basic properties for each point of vortex
        velocity
        gamma
        timePositionX #these 6 will be vectors for graphing of paths over time
        timePositionY
        timePositionZ
        timeVelocityX
        timeVelocityY
        timeVelocityZ
end

#Function to create a new vortex
function newvortex(positionVect, velocityVect, gammaVect)#constructor for type 'vortex.' Used to create a new vortex object
        Vortex(positionVect,velocityVect, gammaVect, [positionVect[1],positionVect[1]], [positionVect[2],positionVect[2]], [positionVect[3],positionVect[3]], [velocityVect[1],velocityVect[1]], [velocityVect[2],velocityVect[2]], [velocityVect[3],velocityVect[3]])
end

#Function that evaluates a stepping method on the vortex objects
function simulatevortex(vortexArray, deltat, endTime)
    t=0.0#time starts here
    nVortecies = length(vortexArray)#VortexArray is an array, each element is a vortex object.
    itteration = 0 #keep track of itterations. start at 0. declare in this scope
    #Main loop for itterating
    while (t < endTime)
        #try maybe doing this without use of round()
        itteration = round(t/deltat + 1)
        #Put current state in vectors
        for n in 1:nVortecies
            #unlike matlab, the existing vectors must be appended with push!
            position =  vortexArray[n].position#get the property "position" from each vortex. Put in an independent 3-element array. This will allow easy access
            push!(vortexArray[n].timePositionX, position[1])#add calculated position to stored chains of data(graphed later)
            push!(vortexArray[n].timePositionY, position[2])
            push!(vortexArray[n].timePositionZ, position[3])

            velocity = vortexArray[n].velocity#same as position, but with velocity
            push!(vortexArray[n].timeVelocityX, velocity[1])#add calculated velocity to stored chains of data(graphed later)
            push!(vortexArray[n].timeVelocityY, velocity[2])
            push!(vortexArray[n].timeVelocityZ, velocity[3])
        end

        #Update positions from velocities
        for n in 1:nVortecies
            # .position holds the position vector for each vortex in vortexArray
            vortexArray[n].position = vortexArray[n].position + (vortexArray[n].velocity * deltat)
        end

        #define velocities based on position
        for n in 1:nVortecies
            #find new velocity for this vortex
            newV = [0.0,0.0,0.0]
            for m in 1:nVortecies #m=other vortex n = current vortex
                if m == n#this means that the vortex being considered is the current one.
                    #do nothing
                else
                    #Define distance between vortecies
                    r = vortexArray[m].position - vortexArray[n].position
                    #Add influence of each vortex
                    newV = newV + (cross(vortexArray[m].gamma, r) ./ (-2.0*pi*dot(r, r))) # +  (G X r)/(2pi r^2) r.r is r^2
                end
            end
            #update velocity for this vortex
            vortexArray[n].velocity = newV
        end

        #update the time
        t = t + deltat


    end #End of while loop

    #Make a time vector to be returned
    timeVector = [-0.001,0.0] #starts out as zero, just like everything else. The "-0.001" just lets the code recognize it as an array.
             #also has itteration + 2 indecies, just like the other graphable vectors
    for i in 1:itteration
        push!(timeVector, (itteration)*deltat)#may need revision, this is why there are parenthesis.
    end

    #Return values for the function
    return (vortexArray,time)
end

#using must be put after function definitions for some reason...
using LinearAlgebra#Used for cross and dot products
#Write the main script here
#initialize an array
testedVortecies = []

#PARAMETERS - read from parameters.simple
using ConfParser
config = ConfParse("parameters.simple")
parse_conf!(config)
#below, the parse function converts output from string to Float64
vortexDiameter = parse(Float64, retrieve(config,"vortexDiameter"))
vortexSeparation = parse(Float64, retrieve(config,"vortexSeparation"))
timeStep = parse(Float64, retrieve(config,"timeStep"))
timeEnd = parse(Float64, retrieve(config,"timeEnd"))




#build the array of vortecies
push!(testedVortecies, newvortex([0.0,0.0,vortexDiameter/2],[0.0,0.0,0.0],[0.0,-1.0,0.0]))
push!(testedVortecies, newvortex([0.0,0.0,-vortexDiameter/2],[0.0,0.0,0.0],[0.0,1.0,0.0]))
push!(testedVortecies, newvortex([vortexSeparation,0.0,vortexDiameter/2],[0.0,0.0,0.0],[0.0,-1.0,0.0]))
push!(testedVortecies, newvortex([vortexSeparation,0.0,-vortexDiameter/2],[0.0,0.0,0.0],[0.0,1.0,0.0]))

(testedVortecies,graphableTime) = simulatevortex(testedVortecies, timeStep, timeEnd)

#necessary to use plots
using Plots

#finds out how many vorticies there are so it knows how many to plot.
plottedVortecies = length(testedVortecies)

#Commented out, used to make a simple plot
#plot3d(testedVortecies[1].timePositionX,testedVortecies[1].timePositionY,testedVortecies[1].timePositionZ)
#plot!(testedVortecies[2].timePositionX,testedVortecies[2].timePositionY,testedVortecies[2].timePositionZ)
#plot!(testedVortecies[3].timePositionX,testedVortecies[3].timePositionY,testedVortecies[3].timePositionZ)
#plot!(testedVortecies[4].timePositionX,testedVortecies[4].timePositionY,testedVortecies[4].timePositionZ)
#multi3Dplot(testedVortecies,plottedVortecies)
#anim = Animation()
#for x = range(0, stop = .5, length = 2)
#    push!(p, x, Float64[sin(x), cos(x)])
#    frame(anim)
#end
animEnd = length(testedVortecies[1].timePositionX)#animate for the same time as number of numerical method itterations were made.

#the code inside this scope could be made better by finding a way to itterativley add plots to a plot3d() with plot! inside a for loop or something.
#Perhaps try this by storing each timePositionX, etc. in an array? see line 39
anim = @animate for i ∈ 1:500:animEnd #start,step,end
    #first ring
    plot3d(testedVortecies[1].timePositionX[1:i],testedVortecies[1].timePositionY[1:i],testedVortecies[1].timePositionZ[1:i], xlims = (0,testedVortecies[1].timePositionX[animEnd] ), ylims = (-vortexDiameter, vortexDiameter), zlims = (-vortexDiameter, vortexDiameter), c = :steelblue)
    plot!(testedVortecies[2].timePositionX[1:i],testedVortecies[2].timePositionY[1:i],testedVortecies[2].timePositionZ[1:i], c = :steelblue)
    #for Ambiance, optional
    #plot!(testedVortecies[1].timePositionX[1:i],testedVortecies[1].timePositionZ[1:i],testedVortecies[1].timePositionY[1:i], c = :steelblue)
    #plot!(testedVortecies[2].timePositionX[1:i],testedVortecies[2].timePositionZ[1:i],testedVortecies[2].timePositionY[1:i], c = :steelblue)
    #circle
    cx = testedVortecies[1].timePositionX[i]#circle position
    cy = testedVortecies[1].timePositionZ[i]#circle radius

    plot!(ones(36).*cx,sin.(0:.2:7).*cy,cos.(0:.2:7).*cy, c = :blue)
    #2nd ring
    plot!(testedVortecies[3].timePositionX[1:i],testedVortecies[3].timePositionY[1:i],testedVortecies[3].timePositionZ[1:i], c = :pink)
    plot!(testedVortecies[4].timePositionX[1:i],testedVortecies[4].timePositionY[1:i],testedVortecies[4].timePositionZ[1:i], c = :pink)
    #for Ambiance, optional
    #plot!(testedVortecies[3].timePositionX[1:i],testedVortecies[3].timePositionZ[1:i],testedVortecies[3].timePositionY[1:i], c = :pink)
    #plot!(testedVortecies[4].timePositionX[1:i],testedVortecies[4].timePositionZ[1:i],testedVortecies[4].timePositionY[1:i], c = :pink)

    #circle
    cx = testedVortecies[3].timePositionX[i]#circle position
    cy = testedVortecies[3].timePositionZ[i]#circle radius

    plot!(ones(36).*cx,sin.(0:.2:7).*cy,cos.(0:.2:7).*cy, c = :red)
end
#save the gif
gif(anim, "Leapfrog.gif", fps = 15)
