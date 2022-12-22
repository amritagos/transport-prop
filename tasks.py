from invoke import task

@task
def build(ctx, install=False):
    print("Building the extension")
    with ctx.cd("transportProp/fastcpp"):
        ctx.run("meson setup bbdir")
        ctx.run("meson compile -C bbdir")
        ctx.run("cp bbdir/*.so .")
    if install:
        ctx.run("flit install")
    else:
        ctx.run("flit build")
