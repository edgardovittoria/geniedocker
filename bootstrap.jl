(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using Geniedocker
const UserApp = Geniedocker
Geniedocker.main()
