from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
# print(config["drobpoxtoken"])
token = open("/home/knut/Projects/dropbox.token").read().strip()
DBox = DropboxRemoteProvider(oauth2_access_token=token)
rule test:
    input:
        DBox.remote("Shoeshine Boys/Real Book Vol 1.pdf")
    shell:
        "echo 'Hello World'"
